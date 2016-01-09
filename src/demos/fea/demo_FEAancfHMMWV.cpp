//
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2013 Project Chrono
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file at the top level of the distribution
// and at http://projectchrono.org/license-chrono.txt.
//

//   Demos code about
//
//     - FEA using ANCF (introduction to dynamics)
//

// Include some headers used by this tutorial...

#include "chrono/physics/ChSystem.h"
#include "chrono/lcp/ChLcpIterativeMINRES.h"
#include "chrono_mkl/ChLcpMklSolver.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono/utils/ChUtilsValidation.h"
#include "chrono/core/ChFileutils.h"
#include "physics/ChLoadContainer.h"
#include "chrono_fea/ChElementShellANCF.h"
#include "chrono_fea/ChElementBrick.h"
#include "chrono_fea/ChMesh.h"
#include "chrono_fea/ChLinkPointFrame.h"
#include "chrono_fea/ChLinkDirFrame.h"
#include "physics/ChLoadContainer.h"  // For Alessandro classes For visualization
#include "chrono_fea/ChVisualizationFEAmesh.h"
#include "chrono/core/ChRealtimeStep.h"
#include "chrono_irrlicht/ChBodySceneNode.h"
#include "chrono_irrlicht/ChBodySceneNodeTools.h"
#include "chrono_irrlicht/ChIrrAppInterface.h"
#include "chrono_irrlicht/ChIrrApp.h"
#include <irrlicht.h>

// Remember to use the namespace 'chrono' because all classes
// of Chrono::Engine belong to this namespace and its children...

using namespace chrono;
using namespace fea;
using namespace irr;
using namespace scene;

class ChCoulombFriction :public ChLoadCustomMultiple{
public:
	ChCoulombFriction(std::vector<ChSharedPtr<ChLoadable>>&mloadables) :ChLoadCustomMultiple(mloadables){};

	int NumContact = 0;
	const double Mu0 = 0.6;
	double ContactLine;
	ChVector<> NetContactForce;


	virtual void ComputeQ(ChState*      state_x, ///< state position to evaluate Q
		ChStateDelta* state_w  ///< state speed to evaluate Q
		) {

		ChMatrixDynamic<double> FlexPos(6, loadables.size());
		ChMatrixDynamic<double> FlexVel(6, loadables.size());
		double Kg;
		double Cg;
		double DeltaDis;
		double DeltaVel;
		double Mu;

		
		if (state_x && state_w) {

			FlexPos.Reset();
			FlexVel.Reset();
			NetContactForce.SetNull();
			for (int ie = 0; ie < loadables.size(); ie++)
			{
				ChVector<> P1 = state_x->ClipVector(6 * ie, 0);
				ChVector<> P1d = state_x->ClipVector(6 * ie + 3, 0);
				ChVector<> V1 = state_w->ClipVector(6 * ie, 0);
				ChVector<> V1d = state_w->ClipVector(6 * ie + 3, 0);

				FlexPos(0, ie) = P1.x;
				FlexPos(1, ie) = P1.y;
				FlexPos(2, ie) = P1.z;
				FlexPos(3, ie) = P1d.x;
				FlexPos(4, ie) = P1d.y;
				FlexPos(5, ie) = P1d.z;

				FlexVel(0, ie) = V1.x;
				FlexVel(1, ie) = V1.y;
				FlexVel(2, ie) = V1.z;
				FlexVel(3, ie) = V1d.x;
				FlexVel(4, ie) = V1d.y;
				FlexVel(5, ie) = V1d.z;

			}

			if (NumContact == 0)
			{
				Kg = 1.0e5;
			}
			else{
				Kg = 8.0e6 / double(NumContact);
			}
			Cg = Kg*0.001;


			for (int ie = 0; ie < loadables.size(); ie++)
			{

				if (FlexPos(2, ie) < ContactLine)
				{
					//GetLog() << i << "\n" << DisFlex(0, i) << "\n" << DisFlex(1, i) << "\n" << DisFlex(2, i) << "\n";

					DeltaDis = FlexPos(2, ie) - ContactLine;
					DeltaVel = FlexVel(2, ie);
					this->load_Q(6 * ie + 2) = -Kg*DeltaDis - Cg*DeltaVel*sqrt(DeltaDis*DeltaDis);


					//Coulomb Friction
					double RelVelX = FlexVel(0, ie);
					double NormVelX = sqrt(RelVelX*RelVelX);
					Mu = Mu0*atan(2.0*RelVelX)*2.0 / CH_C_PI;
					this->load_Q(6 * ie) = -(Mu*(this->load_Q(6 * ie + 2))*RelVelX / (NormVelX));
					//GetLog() << "RelVelX: " << RelVelX << "\nMu: " << Mu << "\nNormVelX: " << NormVelX << "\nFx: " << ContactFRC(0, i) << "\nFz: " << ContactFRC(2, i) << "\natan: " << atan(2.0*RelVelX) << "\n\n";
					if (NormVelX == 0){ this->load_Q(6 * ie) = 0.0; }

					double RelVelY = FlexVel(1, ie);
					double NormVelY = sqrt(RelVelY*RelVelY);
					Mu = Mu0*atan(2.0*RelVelY)*2.0 / CH_C_PI;
					this->load_Q(6 * ie + 1) = -(Mu*(this->load_Q(6 * ie + 2))*RelVelY / (NormVelY));
					if (NormVelY == 0){ this->load_Q(6 * ie + 1) = 0.0; }

					if (this->load_Q(6*ie+2) != 0.0)
					{
						NetContactForce.x += this->load_Q(6 * ie);
						NetContactForce.y += this->load_Q(6 * ie + 1);
						NetContactForce.z += this->load_Q(6 * ie + 2);
					}
				}
			}
		} // end of if(state) loop

	} // end of Compute_Q

	virtual bool IsStiff() { return true; }

};


void ReadInputFile(ChMatrixNM<double, 3000,6> &COORDFlex, ChMatrixNM<double, 3000,6> &VELCYFlex, ChMatrixNM<int, 2880,4> &NodesPerElement, int &TotalNumElements, int &NumElements_x, int &NumElements_y, int &TotalNumNodes, ChMatrixNM<int, 2880,1> &SectionID, ChMatrixNM<double, 10,3> &LayerPROP, ChMatrixNM<int, 10,7> &MatID, ChMatrixNM<double, 7,12> &MPROP, ChMatrixNM<double, 2880,2> &ElementLength)
{
	FILE *inputfile;
	char str1[100];
	int numFlexBody = 0;
	int dummy;
	int count;

	//double ACCELFlex[4000][6];
	//double ACCELRigid[2][7];
	//double LuGreZStart[25][40];
	//double LuGreZStart_dt[25][40];
	//double LuGreZStart_dtdt[25][40];
	double LayPROP[10][7][2];
	int NDR[4000][6];
	int NumLayer[10];
	int MaxSectionNumber = 0;
	int MaxMatID = 0;
	int MTYPE = 0;

	int MAXCOUNT = 100;
	inputfile = fopen("IndataBiLinearShell_Tire(HMMWV).DAT", "r");
	printf("Open IndataBiLinearShell_Tire1.DAT\n");
	if (inputfile == NULL){
		printf("Input data file not found!!\n");
		system("pause");
		exit(1);
	}
	
	TotalNumElements = 0;
	NumElements_x = 0;
	NumElements_y = 0;
	TotalNumNodes = 0;

	//!--------------------------------------!
	//!-- Elememt data            -----------!
	//!--------------------------------------!

	fgets(str1, MAXCOUNT, inputfile);
	printf("%s\n", str1);
	fscanf(inputfile, "%d\n", &numFlexBody);

	fgets(str1, MAXCOUNT, inputfile);
	printf("%s\n", str1);
	fscanf(inputfile, "%d %d %d %d\n", &TotalNumElements, &NumElements_x, &NumElements_y, &TotalNumNodes);
	fgets(str1, MAXCOUNT, inputfile);

	printf("%s\n", str1);
	for (int i = 0; i<TotalNumElements; i++)
	{
		fscanf(inputfile, "%d %d %d %d %d %d %d\n", &count, &dummy, &SectionID(i,0), &NodesPerElement(i, 0), &NodesPerElement(i, 1), &NodesPerElement(i, 2), &NodesPerElement(i, 3));
		printf("SectionID[i] %d\n  ", SectionID(i,0));

		fscanf(inputfile, " %lf %lf\n", &ElementLength(i, 0), &ElementLength(i, 1));
		if (MaxSectionNumber<SectionID(i,0))
		{
			MaxSectionNumber = SectionID(i,0);
		}
		//if(TotalNumNodes<max(NumNodes[i][0],max(NumNodes[i][1],max(NumNodes[i][2],NumNodes[i][3]))))
		//{TotalNumNodes=max(NumNodes[i][0],max(NumNodes[i][1],max(NumNodes[i][2],NumNodes[i][3])));}

		//printf("MaxSectionNumber %lf, %lf \n ", ElemLengthXY[i][0],ElemLengthXY[i][1]);
	}

	//!--------------------------------------!
	//!-- NDR,COORDFlex,VELCYFlex -----------!
	//!--------------------------------------!
	//fscanf(inputfile,"%s\n",str1);
	fgets(str1, MAXCOUNT, inputfile);
	printf("%s\n", str1);
	for (int i = 0; i<TotalNumNodes; i++)
	{
		fscanf(inputfile, "%d %d %d %d %d %d %d\n", &count, &NDR[i][0], &NDR[i][1], &NDR[i][2], &NDR[i][3], &NDR[i][4], &NDR[i][5]);
		fscanf(inputfile, "%lf %lf %lf %lf %lf %lf\n", &COORDFlex(i, 0), &COORDFlex(i, 1), &COORDFlex(i, 2), &COORDFlex(i, 3), &COORDFlex(i, 4), &COORDFlex(i, 5));
		fscanf(inputfile, "%lf %lf %lf %lf %lf %lf\n", &VELCYFlex(i, 0), &VELCYFlex(i, 1), &VELCYFlex(i, 2), &VELCYFlex(i, 3), &VELCYFlex(i, 4), &VELCYFlex(i, 5));
		//printf("NumNodes %d %d %d %d %d %d\n",NDR[i][0],NDR[i][1],NDR[i][2],NDR[i][3],NDR[i][4],NDR[i][5]);
		//printf("NumNodes %lf %lf %lf %lf %lf %lf\n",COORDFlex[i][0],COORDFlex[i][1],COORDFlex[i][2],COORDFlex[i][3],COORDFlex[i][4],COORDFlex[i][5]);
		//system("pause");
	}

	//!--------------------------------------!
	//!--- Read Layer Data ------------------!
	//!--------------------------------------!
	//fscanf(inputfile,"%s\n",str1);
	fgets(str1, MAXCOUNT, inputfile);
	printf("%s\n", str1);
	for (int i = 0; i<MaxSectionNumber; i++)
	{
		fscanf(inputfile, "%d %d\n", &count, &NumLayer[i]);
		for (int j = 0; j<NumLayer[i]; j++)
		{
			fscanf(inputfile, "%lf %lf %d\n", &LayerPROP(MaxSectionNumber*i + j, 1), &LayerPROP(MaxSectionNumber*i + j, 2), &MatID(i, j));
			if (MaxMatID<MatID(i, j))
			{
				MaxMatID = MatID(i, j);
			}
			LayerPROP(MaxSectionNumber*i + j, 0) = NumLayer[i];
			//printf("%lf %lf %d\n",LayPROP[i][j][0],LayPROP[i][j][1],MatID[i][j]);
		}
		//system("pause");
	}

	//!--------------------------------------!
	//!--- Read Material Data ---------------!
	//!--------------------------------------!
	//fscanf(inputfile,"%s\n",str1);
	fgets(str1, MAXCOUNT, inputfile);
	printf("%s\n", str1);
	for (int i = 0; i<MaxMatID; i++)
	{
		fscanf(inputfile, "%d %d\n", &count, &MTYPE);
		if (MTYPE == 1)
		{
			fscanf(inputfile, "%lf %lf %lf %lf\n", &MPROP(i, 0), &MPROP(i, 1), &MPROP(i, 2), &MPROP(i, 3));
		}
		if (MTYPE == 2)
		{
			fscanf(inputfile, "%lf %lf %lf %lf\n", &MPROP(i, 0), &MPROP(i, 1), &MPROP(i, 2), &MPROP(i, 3));
			fscanf(inputfile, "%lf %lf %lf %lf %lf %lf\n", &MPROP(i, 4), &MPROP(i, 5), &MPROP(i, 6), &MPROP(i, 7), &MPROP(i, 8), &MPROP(i, 9));
			
		}
		//printf("%lf %lf %lf %lf\n",MPROP[i][0],MPROP[i][1],MPROP[i][2],MPROP[i][3]);
	}
	
};

void ReadRestartInput()
{
//	FILE *inputfile1;
//	inputfile1 = fopen("QSOL0_All0.txt", "r");
//	printf("Open QSol0.txt\n");
//	if (inputfile1 == NULL){
//		printf("Restart Information not found!\n");
//		system("pause");
//		exit(1);
//	}
//
//	fgets(str1, MAXCOUNT, inputfile1);
//	printf("%s\n", str1);
//	fscanf(inputfile1, "%lf %lf %lf %lf %lf %lf %lf", &COORDRigid[0][0], &COORDRigid[0][1], &COORDRigid[0][2], &COORDRigid[0][3], &COORDRigid[0][4], &COORDRigid[0][5], &COORDRigid[0][6]);
//	fscanf(inputfile1, "%lf %lf %lf %lf %lf %lf %lf", &COORDRigid[1][0], &COORDRigid[1][1], &COORDRigid[1][2], &COORDRigid[1][3], &COORDRigid[1][4], &COORDRigid[1][5], &COORDRigid[1][6]);
//	for (int i = 0; i < TotalNumNodes; i++)
//	{
//		fscanf(inputfile1, "%lf %lf %lf %lf %lf %lf", &COORDFlex[i][0], &COORDFlex[i][1], &COORDFlex[i][2], &COORDFlex[i][3], &COORDFlex[i][4], &COORDFlex[i][5]);
//	}
//	for (int i = 0; i < 25; i++)
//	{
//		for (int j = 0; j < 40; j++)
//		{
//			fscanf(inputfile1, "%lf ", &LuGreZStart[i][j]);
//		}
//	}
//	fscanf(inputfile1, "\n");
//	fgets(str1, MAXCOUNT, inputfile1);
//	printf("%s\n", str1);
//	fscanf(inputfile1, "%lf %lf %lf %lf %lf %lf %lf", &VELCYRigid[0][0], &VELCYRigid[0][1], &VELCYRigid[0][2], &VELCYRigid[0][3], &VELCYRigid[0][4], &VELCYRigid[0][5], &VELCYRigid[0][6]);
//	fscanf(inputfile1, "%lf %lf %lf %lf %lf %lf %lf", &VELCYRigid[1][0], &VELCYRigid[1][1], &VELCYRigid[1][2], &VELCYRigid[1][3], &VELCYRigid[1][4], &VELCYRigid[1][5], &VELCYRigid[1][6]);
//	for (int i = 0; i < TotalNumNodes; i++)
//	{
//		fscanf(inputfile1, "%lf %lf %lf %lf %lf %lf", &VELCYFlex[i][0], &VELCYFlex[i][1], &VELCYFlex[i][2], &VELCYFlex[i][3], &VELCYFlex[i][4], &VELCYFlex[i][5]);
//	}
//	for (int i = 0; i < 25; i++)
//	{
//		for (int j = 0; j < 40; j++)
//		{
//			fscanf(inputfile1, "%lf ", &LuGreZStart_dt[i][j]);
//		}
//	}
//	fscanf(inputfile1, "\n");
//	fgets(str1, MAXCOUNT, inputfile1);
//	printf("%s\n", str1);
//	fscanf(inputfile1, "%lf %lf %lf %lf %lf %lf %lf", &ACCELRigid[0][0], &ACCELRigid[0][1], &ACCELRigid[0][2], &ACCELRigid[0][3], &ACCELRigid[0][4], &ACCELRigid[0][5], &ACCELRigid[0][6]);
//	fscanf(inputfile1, "%lf %lf %lf %lf %lf %lf %lf", &ACCELRigid[1][0], &ACCELRigid[1][1], &ACCELRigid[1][2], &ACCELRigid[1][3], &ACCELRigid[1][4], &ACCELRigid[1][5], &ACCELRigid[1][6]);
//	for (int i = 0; i < TotalNumNodes; i++)
//	{
//		fscanf(inputfile1, "%lf %lf %lf %lf %lf %lf", &ACCELFlex[i][0], &ACCELFlex[i][1], &ACCELFlex[i][2], &ACCELFlex[i][3], &ACCELFlex[i][4], &ACCELFlex[i][5]);
//	}
//	for (int i = 0; i < 25; i++)
//	{
//		for (int j = 0; j < 40; j++)
//		{
//			fscanf(inputfile1, "%lf ", &LuGreZStart_dtdt[i][j]);
//		}
//	}
//	GetLog() << "Restart Complete!\n\n";
};


int main(int argc, char* argv[]) {

	GetLog() << "\n-------------------------------------------------\n";
	GetLog() << "TEST: ANCF Tire (Fixed),  implicit integration \n\n";

	FILE *inputfile; // input stream
	FILE *outputfile;
	FILE *outputfile1;
	FILE *outputfile2;
	FILE *outputfile3;
	FILE *outputfile4;

	// The physical system: it contains all physical objects.
	ChSystem my_system;

	// Create a mesh, that is a container for groups
	// of elements and their referenced nodes.
	ChSharedPtr<ChMesh> my_mesh(new ChMesh);
	
	int TotalNumNodes;
	ChMatrixNM<double, 3000, 6> COORDFlex;
	ChMatrixNM<double, 3000, 6> VELCYFlex;
	ChMatrixNM<double, 3000, 6> ACCELFlex;
	ChMatrixNM<double, 7, 1> COORDRigid;
	ChMatrixNM<double, 7, 1> VELCYRigid;
	ChMatrixNM<double, 7, 1> ACCELRigid;
	ChMatrixNM<int, 2880, 4> NodesPerElement;
	ChMatrixNM<double, 2880, 2> ElemLength;
	int TotalNumElements;
	int NumElements_x;
	int NumElements_y;
	ChMatrixNM<int, 2880,1> SectionID;
	ChMatrixNM<double, 10, 3> LayPROP;
	ChMatrixNM<int, 10, 7> MatID;
	ChMatrixNM<double, 7, 12> MPROP;

	bool Restart = false;


	ReadInputFile(COORDFlex,VELCYFlex,NodesPerElement, TotalNumElements, NumElements_x, NumElements_y,TotalNumNodes, SectionID, LayPROP, MatID, MPROP, ElemLength);


	///////////////////////////////////////////////////////////////////////////
	//// Material List (for HMMWV)
	//// i=0: Carcass
	//// i=1: Steel belt in rubber matrix
	//// i=2: Rubber
	///////////////////////////////////////////////////////////////////////////
	std::vector<ChSharedPtr<ChMaterialShellANCF>> MaterialList(MPROP.GetRows());
		
	for (int i = 0; i < MPROP.GetRows(); i++)
	{
		double rho = MPROP(i, 0);
		ChVector<double> E(MPROP(i, 1), MPROP(i, 2), MPROP(i, 3));
		ChVector<double> nu(MPROP(i, 4), MPROP(i, 5), MPROP(i, 6));
		ChVector<double> G(MPROP(i, 7), MPROP(i, 8), MPROP(i, 9));
		MaterialList[i] = ChSharedPtr<ChMaterialShellANCF>(new ChMaterialShellANCF(rho, E, nu, G));
	}


	for (int i = 0; i < TotalNumNodes; i++)
	{
		ChSharedPtr<ChNodeFEAxyzD> node(new ChNodeFEAxyzD(ChVector<>(COORDFlex(i, 0), COORDFlex(i, 1), COORDFlex(i, 2)), ChVector<>(COORDFlex(i, 3), COORDFlex(i, 4), COORDFlex(i, 5))));
		node->SetPos_dt(ChVector<>(VELCYFlex(i, 0), VELCYFlex(i, 1), VELCYFlex(i, 2)));
		node->SetD_dt(ChVector<>(VELCYFlex(i, 3), VELCYFlex(i, 4), VELCYFlex(i, 5)));
		node->SetPos_dtdt(ChVector<>(ACCELFlex(i, 0), ACCELFlex(i, 1), ACCELFlex(i, 2)));
		node->SetD_dtdt(ChVector<>(ACCELFlex(i, 3), ACCELFlex(i, 4), ACCELFlex(i, 5)));
		node->SetMass(0.0);
		my_mesh->AddNode(node);
	}
	GetLog() << "TotalNumNodes: " << TotalNumNodes << "\n\n";
	ChSharedPtr<ChNodeFEAxyzD> nodetip(my_mesh->GetNode((TotalNumElements / 2)).DynamicCastTo<ChNodeFEAxyzD>());
	GetLog() << "X : " << nodetip->GetPos().x << " Y : " << nodetip->GetPos().y << " Z : " << nodetip->GetPos().z << "\n\n";
	GetLog() << "dX : " << nodetip->GetD().x << " dY : " << nodetip->GetD().y << " dZ : " << nodetip->GetD().z << "\n\n";
	double ContactZ = -0.02;
	int elemcount = 0;
	double timestep = 0.00005;

	for (int i = 0; i < TotalNumElements; i++)
	{
		ChSharedPtr<ChElementShellANCF> element(new ChElementShellANCF);
		element->SetNodes(my_mesh->GetNode(NodesPerElement(i,0) - 1).DynamicCastTo<ChNodeFEAxyzD>(), my_mesh->GetNode(NodesPerElement(i,1) - 1).DynamicCastTo<ChNodeFEAxyzD>(), my_mesh->GetNode(NodesPerElement(i,2) - 1).DynamicCastTo<ChNodeFEAxyzD>(), my_mesh->GetNode(NodesPerElement(i,3) - 1).DynamicCastTo<ChNodeFEAxyzD>());
		element->SetDimensions(ElemLength(i, 0), ElemLength(i, 1));

		for (int j = 0; j < LayPROP(3*(SectionID(i) - 1), 0); j++)
		{
			element->AddLayer(LayPROP(3 * (SectionID(i) - 1) + j, 1), LayPROP(3 * (SectionID(i) - 1) + j, 2) * CH_C_DEG_TO_RAD, MaterialList[MatID(SectionID(i) - 1, j) - 1]);
			//GetLog() << "Thickness: " << LayPROP(3 * (SectionID(i) - 1) + j, 1) << "  Ply: " << LayPROP(3 * (SectionID(i) - 1) + j, 2) << "  Mat: " << MatID(SectionID(i) - 1, j)<<"\n";
		}
		element->SetAlphaDamp(0.01);
		element->SetGravityOn(true);
		my_mesh->AddElement(element);
	}

	ChSharedPtr<ChBody> Rim(new ChBody);
	my_system.Add(Rim);
	Rim->SetBodyFixed(false);
	Rim->SetPos(ChVector<>(0.0, 0.0, 0.4673)); //0.46736
	Rim->SetRot(ChQuaternion<>(1.0, 0.0, 0.0, 0.0));
	Rim->SetPos_dt(ChVector<>(0.0, 0.0, 0.0));
	Rim->SetRot_dt(ChQuaternion<>(0.0, 0.0, 0.0, 0.0));
	Rim->SetMass(509.69);
	Rim->SetInertiaXX(ChVector<>(0.1457270, 0.2359220, 0.1457270));
	if (Restart)
	{
		Rim->SetPos(ChVector<>(COORDRigid(0, 0), COORDRigid(1, 0), COORDRigid(2, 0)));
		Rim->SetRot(ChQuaternion<>(COORDRigid(3, 0), COORDRigid(4, 0), COORDRigid(5, 0), COORDRigid(6, 0)));
		Rim->SetPos_dt(ChVector<>(VELCYRigid(0, 0), VELCYRigid(1, 0), VELCYRigid(2, 0)));
		Rim->SetRot_dt(ChQuaternion<>(VELCYRigid(3, 0), VELCYRigid(4, 0), VELCYRigid(5, 0), VELCYRigid(5, 0)));
		Rim->SetPos_dtdt(ChVector<>(ACCELRigid(0, 0), ACCELRigid(1, 0), ACCELRigid(2, 0)));
		Rim->SetRot_dtdt(ChQuaternion<>(ACCELRigid(3, 0), ACCELRigid(4, 0), ACCELRigid(5, 0), ACCELRigid(6, 0)));
	}

	ChSharedPtr<ChBody> Ground(new ChBody);
	Ground->SetBodyFixed(true);
	Ground->SetPos(ChVector<>(0.0, 0.0, -0.02));
	Ground->SetRot(ChQuaternion<>(1.0, 0.0, 0.0, 0.0));
	my_system.Add(Ground);

	my_system.Set_G_acc(ChVector<>(0, 0, -9.81));

	ChSharedPtr<ChLinkPointFrame> constraint;
	ChSharedPtr<ChLinkDirFrame> constraintD;
	ChSharedPtr<ChNodeFEAxyzD> ConstrainedNode;
	ChSharedPtr<ChLinkLockPlanePlane> constraintRim;

	// First: loads must be added to "load containers",
	// and load containers must be added to your ChSystem
	ChSharedPtr<ChLoadContainer> Mloadcontainer(new ChLoadContainer);

	//// LuGre Load Class initialization
	std::vector<ChSharedPtr<ChCoulombFriction>> CFLoadList(NumElements_y+1);
	for (int i = 0; i < NumElements_y + 1; i++)
	{

		std::vector< ChSharedPtr< ChLoadable > > mnodelist;
		for (int inode = 0; inode < (TotalNumNodes / (NumElements_y + 1)); inode++)
		{
			ChSharedPtr<ChNodeFEAxyzD> FrictionNode(my_mesh->GetNode((i*TotalNumNodes / (NumElements_y + 1)) + inode).DynamicCastTo<ChNodeFEAxyzD>());
			mnodelist.push_back(FrictionNode);
		}

		CFLoadList[i] = ChSharedPtr<ChCoulombFriction>(new ChCoulombFriction(mnodelist));
		CFLoadList[i]->ContactLine = ContactZ;
		Mloadcontainer->Add(CFLoadList[i]);
	}

	class MyPressureLoad : public ChLoaderUVdistributed {
	private:
		ChSharedPtr<ChElementShellANCF> m_element;

	public:
		// Useful: a constructor that also sets ChLoadable
		MyPressureLoad(ChSharedPtr<ChLoadableUV> element) : ChLoaderUVdistributed(element) {
			m_element = element.StaticCastTo<ChElementShellANCF>();
		};
		virtual bool IsStiff() override { return true; }
		virtual int GetIntegrationPointsU() { return 2; }
		virtual int GetIntegrationPointsV() { return 2; }
		// Compute F=F(u)
		// This is the function that you have to implement. It should return the
		// load at U. For Eulero beams, loads are expected as 6-rows vectors, containing
		// a wrench: forceX, forceY, forceZ, torqueX, torqueY, torqueZ.
		ChVector<> FPressure;
		virtual void ComputeF(
			const double U,
			const double V,              ///< parametric coordinate in line
			ChVectorDynamic<>& F,        ///< Result F vector here, size must be = n.field coords.of loadable
			ChVectorDynamic<>* state_x,  ///< if != 0, update state (pos. part) to this, then evaluate F
			ChVectorDynamic<>* state_w   ///< if != 0, update state (speed part) to this, then evaluate F
			) {
			ChVector<> Position1;
			ChVector<> Gradient1;
			ChVector<> Position2;
			ChVector<> Gradient2;
			ChVector<> Position3;
			ChVector<> Gradient3;
			ChVector<> Position4;
			ChVector<> Gradient4;
			double PressureVal = 220e3;  // Pressure

			if (state_x && state_w) {
				Position1 = state_x->ClipVector(0, 0);
				Gradient1 = state_x->ClipVector(3, 0);
				Position2 = state_x->ClipVector(6, 0);
				Gradient2 = state_x->ClipVector(9, 0);
				Position3 = state_x->ClipVector(12, 0);
				Gradient3 = state_x->ClipVector(15, 0);
				Position4 = state_x->ClipVector(18, 0);
				Gradient4 = state_x->ClipVector(21, 0);

				ChMatrixNM<double, 1, 8> Nx;
				ChMatrixNM<double, 1, 8> N;
				ChMatrixNM<double, 1, 8> Ny;
				ChMatrixNM<double, 1, 8> Nz;
				ChMatrixNM<double, 3, 8> d;
				(d).PasteVector(Position1, 0, 0);
				(d).PasteVector(Gradient1, 0, 1);
				(d).PasteVector(Position2, 0, 2);
				(d).PasteVector(Gradient2, 0, 3);
				(d).PasteVector(Position3, 0, 4);
				(d).PasteVector(Gradient3, 0, 5);
				(d).PasteVector(Position4, 0, 6);
				(d).PasteVector(Gradient4, 0, 7);
				m_element->ShapeFunctions(N, U, V, 0);
				m_element->ShapeFunctionsDerivativeX(Nx, U, V, 0);
				m_element->ShapeFunctionsDerivativeY(Ny, U, V, 0);
				m_element->ShapeFunctionsDerivativeZ(Nz, U, V, 0);

				ChMatrixNM<double, 1, 3> Nx_d;
				Nx_d.MatrMultiplyT(Nx, d);

				ChMatrixNM<double, 1, 3> Ny_d;
				Ny_d.MatrMultiplyT(Ny, d);

				ChMatrixNM<double, 1, 3> Nz_d;
				Nz_d.MatrMultiplyT(Nz, d);

				ChMatrixNM<double, 3, 3> rd;
				rd(0, 0) = Nx_d(0, 0);
				rd(1, 0) = Nx_d(0, 1);
				rd(2, 0) = Nx_d(0, 2);
				rd(0, 1) = Ny_d(0, 0);
				rd(1, 1) = Ny_d(0, 1);
				rd(2, 1) = Ny_d(0, 2);
				rd(0, 2) = Nz_d(0, 0);
				rd(1, 2) = Nz_d(0, 1);
				rd(2, 2) = Nz_d(0, 2);

				ChVector<> G1xG2;
				G1xG2(0) = rd(1, 0) * rd(2, 1) - rd(2, 0) * rd(1, 1);
				G1xG2(1) = rd(2, 0) * rd(0, 1) - rd(0, 0) * rd(2, 1);
				G1xG2(2) = rd(0, 0) * rd(1, 1) - rd(1, 0) * rd(0, 1);
				G1xG2.Normalize();
				FPressure = G1xG2 * PressureVal;
			}
			F.PasteVector(FPressure, 0, 0);
		}
	};

	// Create the load (and handle it with a shared pointer).
	// The ChLoad is a 'container' for your ChLoader.
	// It is created using templates, that is instancing a ChLoad<a_loader_class>()
	// initiate for loop for all the elements
	for (int NoElmPre = 0; NoElmPre < TotalNumElements; NoElmPre++) {
		ChSharedPtr<ChLoad<MyPressureLoad>> PressureElement(
			new ChLoad<MyPressureLoad>(my_mesh->GetElement(NoElmPre).StaticCastTo<ChElementShellANCF>()));
		Mloadcontainer->Add(PressureElement);  // do not forget to add the load to the load container.
	}

	my_system.Add(Mloadcontainer);

	// Remember to add the mesh to the system!
	my_system.Add(my_mesh);

	// Constrain the flexible tire to the rigid Rim body.
	for (int i = 0; i < TotalNumNodes;i++)
	{
		if (i<NumElements_x || i >= TotalNumNodes - NumElements_x){  // Only constrain the nodes at the ends of the bead section

			ConstrainedNode = ChSharedPtr<ChNodeFEAxyzD>(my_mesh->GetNode(i).DynamicCastTo<ChNodeFEAxyzD>());

			// Add position constraints
			constraint = ChSharedPtr<ChLinkPointFrame>(new ChLinkPointFrame);
			constraint->Initialize(ConstrainedNode, Rim);
			my_system.Add(constraint);

			// Add rotation constraints
			constraintD = ChSharedPtr<ChLinkDirFrame>(new ChLinkDirFrame);
			constraintD->Initialize(ConstrainedNode, Rim);
			constraintD->SetDirectionInAbsoluteCoords(ConstrainedNode->GetD());
			my_system.Add(constraintD);
		}
	}
	// Constrain the Rim to the X-Z plane
	constraintRim = ChSharedPtr<ChLinkLockPlanePlane>(new ChLinkLockPlanePlane);
	my_system.AddLink(constraintRim);
	constraintRim->Initialize(Rim, Ground, ChCoordsys<>(ChVector<>(0.0, 0.0, 0.0), Q_from_AngX(CH_C_PI_2)));


	// This is mandatory
	my_system.SetupInitial();

	// Use the MKL Solver
	ChLcpMklSolver* mkl_solver_stab = new ChLcpMklSolver;
	ChLcpMklSolver* mkl_solver_speed = new ChLcpMklSolver;
	my_system.ChangeLcpSolverStab(mkl_solver_stab);
	my_system.ChangeLcpSolverSpeed(mkl_solver_speed);
	mkl_solver_stab->SetSparsityPatternLock(true);
	mkl_solver_speed->SetSparsityPatternLock(true);
	my_system.Update();

	// Choose the time integrator
	my_system.SetIntegrationType(ChSystem::INT_HHT);  
	ChSharedPtr<ChTimestepperHHT> mystepper = my_system.GetTimestepper().DynamicCastTo<ChTimestepperHHT>();
	mystepper->SetAlpha(-0.2);
	mystepper->SetMaxiters(20);
	mystepper->SetTolerance(5e-5);
	mystepper->SetMode(ChTimestepperHHT::POSITION);
	mystepper->SetVerbose(true);
	mystepper->SetScaling(true);

	my_system.Setup();
	my_system.Update();

	// Create output files
	outputfile = fopen("OutPutBody1.txt", "w");				// Time history of nodal coordinates
	outputfile1 = fopen("AllPositions.txt", "w");			// Time history of rigid bodies
	outputfile2 = fopen("ContactInformation.txt", "w");		// Ground line and normal contact force
	outputfile3 = fopen("RES.txt", "w");					// Number of iterations and CPU time
	outputfile4 = fopen("StockContactForce.txt", "w");		// Time history of contact forces per nodal coordinate

	GetLog() << "Contact Line: " << ContactZ << "\n";
	ChVector<> NetContact;
	int NumContact1 = 0;

	GetLog() << "Contact ForceX: " << NetContact.x << "\n";
	GetLog() << "Contact ForceY: " << NetContact.y << "\n";
	GetLog() << "Contact ForceZ: " << NetContact.z << "\n";
	GetLog() << "Number of Contact Points: " << NumContact1 << "\n";
	GetLog() << " t=  " << my_system.GetChTime() << "\n\n";
	fprintf(outputfile, "%15.7e  ", my_system.GetChTime());
	for (int ii = 0; ii < TotalNumNodes; ii++)
	{
		ChSharedPtr<ChNodeFEAxyzD> nodetip(my_mesh->GetNode((ii)).DynamicCastTo<ChNodeFEAxyzD>());
		fprintf(outputfile, "%15.7e  ", nodetip->GetPos().x);
		fprintf(outputfile, "%15.7e  ", nodetip->GetPos().y);
		fprintf(outputfile, "%15.7e  ", nodetip->GetPos().z);
		fprintf(outputfile, "%15.7e  ", nodetip->GetD().x);
		fprintf(outputfile, "%15.7e  ", nodetip->GetD().y);
		fprintf(outputfile, "%15.7e  ", nodetip->GetD().z);
	}
	fprintf(outputfile4, "%15.7e  ", my_system.GetChTime());
	for (int ii = 0; ii < TotalNumNodes; ii++)
	{
		ChSharedPtr<ChNodeFEAxyzD> nodetip(my_mesh->GetNode((ii)).DynamicCastTo<ChNodeFEAxyzD>());
		fprintf(outputfile, "%15.7e  ", nodetip->GetPos_dt().x);
		fprintf(outputfile, "%15.7e  ", nodetip->GetPos_dt().y);
		fprintf(outputfile, "%15.7e  ", nodetip->GetPos_dt().z);
		fprintf(outputfile, "%15.7e  ", nodetip->GetD_dt().x);
		fprintf(outputfile, "%15.7e  ", nodetip->GetD_dt().y);
		fprintf(outputfile, "%15.7e  ", nodetip->GetD_dt().z);
	}

	fprintf(outputfile, "\n  ");
	fprintf(outputfile4, "\n  ");
	fprintf(outputfile1, "%15.7e  ", my_system.GetChTime());
	fprintf(outputfile1, "%15.7e  ", 0.0);
	fprintf(outputfile1, "%15.7e  ", 0.0);
	fprintf(outputfile1, "%15.7e  ", 0.0);
	fprintf(outputfile1, "%15.7e  ", 1.0);
	fprintf(outputfile1, "%15.7e  ", 0.0);
	fprintf(outputfile1, "%15.7e  ", 0.0);
	fprintf(outputfile1, "%15.7e  ", 0.0);
	fprintf(outputfile1, "%15.7e  ", Rim->GetPos().x);
	fprintf(outputfile1, "%15.7e  ", Rim->GetPos().y);
	fprintf(outputfile1, "%15.7e  ", Rim->GetPos().z);
	fprintf(outputfile1, "%15.7e  ", Rim->GetRot().e0);
	fprintf(outputfile1, "%15.7e  ", Rim->GetRot().e1);
	fprintf(outputfile1, "%15.7e  ", Rim->GetRot().e2);
	fprintf(outputfile1, "%15.7e  ", Rim->GetRot().e3);
	fprintf(outputfile1, "\n  ");
	int NumCont = 0;
	//////////////////////////////////////////////////////
	double start = std::clock();
	while (my_system.GetChTime() < 0.2) {

		Rim->Empty_forces_accumulators();
		//Rim->Accumulate_force(ChVector<>(0.0, 0.0, -5000.0), Rim->GetPos(), 0);
		//Rim->Set_Scr_force(ChVector<>(0.0, 0.0, -4000.0));
		//if (/*my_system.GetChTime()>0.09&&*/my_system.GetChTime()+1.9 <= 0.3&&Rim->GetPos_dt().x > -22.571)
		//{
		// Rim->Accumulate_torque(ChVector<>(0.0, -1000.0*(my_system.GetChTime()+1.9), 0.0), 1);
		// //Rim->Set_Scr_torque(ChVector<>(0.0, 100.0*(my_system.GetChTime()), 0.0));
		//}
		//else if (my_system.GetChTime()+1.9>0.3&&my_system.GetChTime()+1.9 <= 0.7&&Rim->GetPos_dt().x > -22.571)
		//{
		// Rim->Accumulate_torque(ChVector<>(0.0, -2000.0*(my_system.GetChTime()+1.9), 0.0), 1);
		// //Rim->Set_Scr_torque(ChVector<>(0.0, 200.0*(my_system.GetChTime()), 0.0));
		//}
		//else if (my_system.GetChTime()+1.9>0.7&&my_system.GetChTime()+1.9 <= 1.5&&Rim->GetPos_dt().x > -22.571)
		//{
		// Rim->Accumulate_torque(ChVector<>(0.0, -3000.0*(my_system.GetChTime()+1.9), 0.0), 1);
		// //Rim->Set_Scr_torque(ChVector<>(0.0, 300.0*(my_system.GetChTime()), 0.0));
		//}
		//else if (my_system.GetChTime()+1.9>1.5&&my_system.GetChTime()+1.9 <= 2.0&&Rim->GetPos_dt().x > -22.571)
		//{
		// Rim->Accumulate_torque(ChVector<>(0.0, -4000.0*(my_system.GetChTime()+1.9), 0.0), 1);
		// //Rim->Set_Scr_torque(ChVector<>(0.0, 400.0*(my_system.GetChTime()), 0.0));
		//}
		//GetLog() << "Acc Force: " << Rim->Get_accumulated_force().x << " " << Rim->Get_accumulated_force().y << " " << Rim->Get_accumulated_force().z << "\nAcc Torque: " << Rim->Get_accumulated_torque().x << " " << Rim->Get_accumulated_torque().y << " " << Rim->Get_accumulated_torque().z << "\n\n";
		
		//==Start analysis==//
		my_system.DoStepDynamics(timestep);

		//==============================//
		//== Output programs ===========//
		//==============================//

		GetLog() << " t=  " << my_system.GetChTime() << "\n\n";
		NumCont = 0;
		NetContact.x = 0.0;
		NetContact.y = 0.0;
		NetContact.z = 0.0;
		fprintf(outputfile, "%15.7e  ", my_system.GetChTime());
		for (int ii = 0; ii < TotalNumNodes; ii++)
		{
			ChSharedPtr<ChNodeFEAxyzD> nodetip(my_mesh->GetNode((ii)).DynamicCastTo<ChNodeFEAxyzD>());
			fprintf(outputfile, "%15.7e  ", nodetip->GetPos().x);
			fprintf(outputfile, "%15.7e  ", nodetip->GetPos().y);
			fprintf(outputfile, "%15.7e  ", nodetip->GetPos().z);
			fprintf(outputfile, "%15.7e  ", nodetip->GetD().x);
			fprintf(outputfile, "%15.7e  ", nodetip->GetD().y);
			fprintf(outputfile, "%15.7e  ", nodetip->GetD().z);
			if (nodetip->GetPos().z < ContactZ)
			{
				NumCont++;
			}
		}
		for (int i = 0; i < NumElements_y + 1; i++)
		{
			CFLoadList[i]->NumContact = NumCont;
			NetContact += CFLoadList[i]->NetContactForce;
		}
		fprintf(outputfile4, "%15.7e  ", my_system.GetChTime());
		for (int ii = 0; ii < TotalNumNodes; ii++)
		{
			ChSharedPtr<ChNodeFEAxyzD> nodetip(my_mesh->GetNode((ii)).DynamicCastTo<ChNodeFEAxyzD>());
			fprintf(outputfile, "%15.7e  ", nodetip->GetPos_dt().x);
			fprintf(outputfile, "%15.7e  ", nodetip->GetPos_dt().y);
			fprintf(outputfile, "%15.7e  ", nodetip->GetPos_dt().z);
			fprintf(outputfile, "%15.7e  ", nodetip->GetD_dt().x);
			fprintf(outputfile, "%15.7e  ", nodetip->GetD_dt().y);
			fprintf(outputfile, "%15.7e  ", nodetip->GetD_dt().z);
		}

		fprintf(outputfile, "\n  ");
		fprintf(outputfile4, "\n   ");
		fprintf(outputfile1, "%15.7e  ", my_system.GetChTime());
		fprintf(outputfile1, "%15.7e  ", 0.0);
		fprintf(outputfile1, "%15.7e  ", 0.0);
		fprintf(outputfile1, "%15.7e  ", 0.0);
		fprintf(outputfile1, "%15.7e  ", 1.0);
		fprintf(outputfile1, "%15.7e  ", 0.0);
		fprintf(outputfile1, "%15.7e  ", 0.0);
		fprintf(outputfile1, "%15.7e  ", 0.0);
		fprintf(outputfile1, "%15.7e  ", Rim->GetPos().x);
		fprintf(outputfile1, "%15.7e  ", Rim->GetPos().y);
		fprintf(outputfile1, "%15.7e  ", Rim->GetPos().z);
		fprintf(outputfile1, "%15.7e  ", Rim->GetRot().e0);
		fprintf(outputfile1, "%15.7e  ", Rim->GetRot().e1);
		fprintf(outputfile1, "%15.7e  ", Rim->GetRot().e2);
		fprintf(outputfile1, "%15.7e  ", Rim->GetRot().e3);
		fprintf(outputfile1, "\n  ");

		fprintf(outputfile2, "%15.7e  ", my_system.GetChTime());
		fprintf(outputfile2, "%15.7e  ", ContactZ);
		fprintf(outputfile2, "%15.7e  ", NetContact.z);
		fprintf(outputfile2, "%d  ", NumCont);
		fprintf(outputfile2, "\n  ");

		fprintf(outputfile3, "%15.7e  ", my_system.GetChTime());
		fprintf(outputfile3, "%d  ", mystepper->GetNumIterations());
		fprintf(outputfile3, "\n  ");

		GetLog() << "Contact Line: " << ContactZ << "\n";
		GetLog() << "Contact ForceX: " << NetContact.x << "\n";
		GetLog() << "Contact ForceY: " << NetContact.y << "\n";
		GetLog() << "Contact ForceZ: " << NetContact.z << "\n";
		GetLog() << "Rim X: " << Rim->GetPos().x << "\n";
		GetLog() << "Rim Z: " << Rim->GetPos().z << "\n";
		GetLog() << "Rim Vel X: " << Rim->GetPos_dt().x << "\n";
		GetLog() << "Number of Contact Points: " << NumCont << "\n\n\n";
	}

	double duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;

	ChSharedPtr<ChTimestepperHHT> mystepper1 = my_system.GetTimestepper().DynamicCastTo<ChTimestepperHHT>();
	GetLog() << "Simulation Time: " << duration << "\n";
	fprintf(outputfile3, "%15.7e  ", duration);
	fprintf(outputfile3, "%d  ", mystepper1->GetNumIterations());
	fprintf(outputfile3, "\n  ");

	system("pause");
	return 0;
}