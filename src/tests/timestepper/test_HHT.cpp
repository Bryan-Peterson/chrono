#include <math.h>
#include <stdio.h>

#include "chrono/core/ChLinearAlgebra.h"
#include "chrono/core/ChLog.h"
#include "chrono/timestepper/ChTimestepper.h"

using namespace chrono;

// 2nd order odcillator problem definition.
// Define a class inherited from ChIntegrableIIorder, which will represent the differential equations
// by implementing the interfaces to implicit solvers.
// We assume   M*a = F(x,v,t)
class OscillatorProblem : public ChIntegrableIIorder {
  private:
    double M;
    double K;
    double R;
    double mT;
    double mx;
    double mv;

  public:
    OscillatorProblem() {
        mT = 0;
        M = 1;
        K = 30;
        R = 0;
        mx = 0;
        mv = 0.6;
    }

    /// the number of coordinates in the state, x position part:
    virtual int GetNcoords_x() { return 1; }

    /// system -> state
    virtual void StateGather(ChState& x, ChStateDelta& v, double& T) {
        x(0) = mx;
        v(0) = mv;
        T = mT;
    }

    /// state -> system
    virtual void StateScatter(const ChState& x, const ChStateDelta& v, const double T) {
        mx = x(0);
        mv = v(0);
        mT = T;
    }

    /// compute  dy/dt=f(y,t)
    virtual void StateSolveA(ChStateDelta& dvdt,              ///< result: computed accel. a=dv/dt
                             ChVectorDynamic<>& L,            ///< result: computed lagrangian multipliers, if any
                             const ChState& x,                ///< current state, x
                             const ChStateDelta& v,           ///< current state, v
                             const double T,                  ///< current time T
                             const double dt,                 ///< timestep (if needed)
                             bool force_state_scatter = true  ///< if false, y and T are not scattered to the
                             /// system, assuming that someone has done
                             /// StateScatter just before
                             ) {
        if (force_state_scatter)
            StateScatter(x, v, T);
        double F = sin(mT * 20) * 0.02;
        dvdt(0) = (1. / M) * (F - K * mx - R * mv);
    }

    /// Compute the correction with linear system
    ///  Dv = [ c_a*M + c_v*dF/dv + c_x*dF/dx ]^-1 * R
    virtual void StateSolveCorrection(
        ChStateDelta& Dv,                ///< result: computed Dv
        ChVectorDynamic<>& L,            ///< result: computed lagrangian multipliers, if any
        const ChVectorDynamic<>& R,      ///< the R residual
        const ChVectorDynamic<>& Qc,     ///< the Qc residual
        const double c_a,                ///< the factor in c_a*M
        const double c_v,                ///< the factor in c_v*dF/dv
        const double c_x,                ///< the factor in c_x*dF/dv
        const ChState& x,                ///< current state, x part
        const ChStateDelta& v,           ///< current state, v part
        const double T,                  ///< current time T
        bool force_state_scatter = true  ///< if false, x,v and T are not scattered to the system, assuming that
        /// someone has done StateScatter just before
        ) {
        if (force_state_scatter)
            this->StateScatter(x, v, T);

        Dv(0) = R(0) * 1.0 / (c_a * this->M + c_v * (-this->R) + c_x * (-this->K));
    }

    ///    R += c*F
    void LoadResidual_F(ChVectorDynamic<>& R,  ///< result: the R residual, R += c*F
                        const double c         ///< a scaling factor
                        ) {
        R(0) += c * (sin(mT * 20) * 0.02 - this->K * mx - this->R * mv);
    }

    ///    R += c*M*w
    void LoadResidual_Mv(ChVectorDynamic<>& R,        ///< result: the R residual, R += c*M*v
                         const ChVectorDynamic<>& w,  ///< the w vector
                         const double c               ///< a scaling factor
                         ) {
        R(0) += c * this->M * w(0);
    }

    /// nothing to do here- no constraints
    virtual void LoadResidual_CqL(ChVectorDynamic<>& R,        ///< result: the R residual, R += c*Cq'*L
                                  const ChVectorDynamic<>& L,  ///< the L vector
                                  const double c               ///< a scaling factor
                                  ) {}

    /// nothing to do here- no constraints
    virtual void LoadConstraint_C(ChVectorDynamic<>& Qc,        ///< result: the Qc residual, Qc += c*C
                                  const double c,               ///< a scaling factor
                                  const bool do_clamp = false,  ///< enable optional clamping of Qc
                                  const double mclam = 1e30     ///< clamping value
                                  ) {}

    /// nothing to do here- no constraints
    virtual void LoadConstraint_Ct(ChVectorDynamic<>& Qc,  ///< result: the Qc residual, Qc += c*Ct
                                   const double c          ///< a scaling factor
                                   ) {}
};

// Pendulum DAE problem.
// Define a class inherited from ChIntegrableIIorder which will represent the differential/algebraic equations
// by implementing the interfaces to implicit solvers.
// We assume   M*a = F(x,v,t)
//             C(x,t)=0;
class PendulumProblem : public ChIntegrableIIorderEasy {
  private:
    double M;
    double K;
    double R;
    double mT;
    double mpx;
    double mpy;
    double mvx;
    double mvy;
    double mlength;
    double mreaction;

  public:
    PendulumProblem() {
        mlength = 1;
        M = 1;
        K = 2;
        R = 0;
        mT = 0;
        mpx = 0;
        mpy = -mlength;
        mvx = 0.2;
        mvy = 0;
        mreaction = 0;
    }

    /// Set pendulum length
    void SetLength(double len) { mlength = len; }

    /// Set (consistent) initial values
    void SetInitialConditions(double x, double y, double vx, double vy) {
        assert(x * x + y * y - mlength * mlength < 1e-10);
        assert(x * vx + y * vy < 1e-10);
        mpx = x;
        mpy = y;
        mvx = vx;
        mvy = vy;
    }

    /// the number of coordinates in the state, x position part:
    virtual int GetNcoords_x() { return 2; }

    /// Tells the number of lagrangian multipliers (constraints)
    virtual int GetNconstr() { return 1; }

    /// system -> state
    virtual void StateGather(ChState& x, ChStateDelta& v, double& T) {
        x(0) = mpx;
        x(1) = mpy;
        v(0) = mvx;
        v(1) = mvy;
        T = mT;
    }

    /// state -> system
    virtual void StateScatter(const ChState& x, const ChStateDelta& v, const double T) {
        mpx = x(0);
        mpy = x(1);
        mvx = v(0);
        mvy = v(1);
        mT = T;
    }

    /// Some timesteppers exploit persistence of reaction information
    virtual void StateGatherReactions(ChVectorDynamic<>& L) { L(0) = mreaction; };
    virtual void StateScatterReactions(const ChVectorDynamic<>& L) { mreaction = L(0); };

    /// Compute the correction with linear system
    ///  Dv = [ c_a*M + c_v*dF/dv + c_x*dF/dx ]^-1 * R
    virtual void StateSolveCorrection(
        ChStateDelta& Dv,                ///< result: computed Dv
        ChVectorDynamic<>& L,            ///< result: computed lagrangian multipliers, if any
        const ChVectorDynamic<>& R,      ///< the R residual
        const ChVectorDynamic<>& Qc,     ///< the Qc residual
        const double c_a,                ///< the factor in c_a*M
        const double c_v,                ///< the factor in c_v*dF/dv
        const double c_x,                ///< the factor in c_x*dF/dv
        const ChState& x,                ///< current state, x part
        const ChStateDelta& v,           ///< current state, v part
        const double T,                  ///< current time T
        bool force_state_scatter = true  ///< if false, x,v and T are not scattered to the system, assuming that
        /// someone has done StateScatter just before
        ) {
        if (force_state_scatter)
            this->StateScatter(x, v, T);

        ChVector<> dirpend(-mpx, -mpy, 0);
        dirpend.Normalize();
        ChVectorDynamic<> b(3);
        b(0) = R(0);
        b(1) = R(1);
        b(2) = Qc(0);
        ChMatrixDynamic<> A(3, 3);
        A(0, 0) = c_a * this->M + c_v * (-this->R) + c_x * (-this->K);
        A(1, 1) = c_a * this->M;
        A(0, 2) = dirpend.x;
        A(1, 2) = dirpend.y;
        A(2, 0) = dirpend.x;
        A(2, 1) = dirpend.y;
        ChVectorDynamic<> w(3);
        ChLinearAlgebra::Solve_LinSys(A, &b, &w);
        Dv(0) = w(0);
        Dv(1) = w(1);
        L(0) = -w(2);  // note assume result sign in multiplier is flipped
    }

    ///    R += c*F
    void LoadResidual_F(ChVectorDynamic<>& R,  ///< result: the R residual, R += c*F
                        const double c         ///< a scaling factor
                        ) {
        R(0) += c * (sin(mT * 20) * 0.000 - this->K * mpx - this->R * mvx);
        R(1) += c * (mT < 0.2 ? (-5) : (50));  // vertical force
    }

    ///    R += c*M*w
    void LoadResidual_Mv(ChVectorDynamic<>& R,        ///< result: the R residual, R += c*M*v
                         const ChVectorDynamic<>& w,  ///< the w vector
                         const double c               ///< a scaling factor
                         ) {
        R(0) += c * this->M * w(0);
        R(1) += c * this->M * w(1);
    }

    ///   R += Cq'*l
    virtual void LoadResidual_CqL(ChVectorDynamic<>& R,        ///< result: the R residual, R += c*Cq'*L
                                  const ChVectorDynamic<>& L,  ///< the L vector
                                  const double c               ///< a scaling factor
                                  ) {
        ChVector<> dirpend(-mpx, -mpy, 0);
        dirpend.Normalize();
        R(0) += c * dirpend.x * L(0);
        R(1) += c * dirpend.y * L(0);
    }

    ///  Qc += c * C
    virtual void LoadConstraint_C(ChVectorDynamic<>& Qc,        ///< result: the Qc residual, Qc += c*C
                                  const double c,               ///< a scaling factor
                                  const bool do_clamp = false,  ///< enable optional clamping of Qc
                                  const double mclam = 1e30     ///< clamping value
                                  ) {
        ChVector<> distpend(-mpx, -mpy, 0);
        Qc(0) += -c * (-distpend.Length() + mlength);
    }

    /// nothing to do here- no rheonomic part
    virtual void LoadConstraint_Ct(ChVectorDynamic<>& Qc,  ///< result: the Qc residual, Qc += c*Ct
                                   const double c          ///< a scaling factor
                                   ) {}
};

// ==========================================================================================================

void Oscillator() {
    printf("\nIntegrate 2nd order oscillator: M*x'' + R*x' + K*x = F(t)\n");

    double timestep = 0.05;
    int num_steps = 100;

    // Create the problem object
    OscillatorProblem mintegrable;

    // Create the HHT timestepper
    ChTimestepperHHT mystepper(&mintegrable);
    mystepper.SetAlpha(0);
    mystepper.SetTolerance(1e-10);
    mystepper.SetMaxiters(6);
    mystepper.SetMaxItersSuccess(3);
    mystepper.SetRequiredSuccessfulSteps(5);
    mystepper.SetStepIncreaseFactor(2);
    mystepper.SetStepDecreaseFactor(0.5);
    mystepper.SetVerbose(true);

    // Execute the time integration
    for (int i = 0; i < num_steps; i++) {
        mystepper.Advance(timestep);
        printf("    %7.4f | %12.8f  %12.8f | %2d\n",
               mystepper.GetTime(),          // time
               mystepper.get_X()(0),         // position
               mystepper.get_V()(0),         // velocity
               mystepper.GetNumIterations()  // number of HHT iterations
               );
    }
}

void Pendulum() {
    printf("\nIntegrate pendulum DAE\n");

    // Solver settings.
    double timestep = 0.05;
    int num_steps = 10;

    // Create the problem object
    PendulumProblem mintegrable;

    // Create the HHT timestepper
    ChTimestepperHHT mystepper(&mintegrable);
    mystepper.SetAlpha(-0.2);
    mystepper.SetTolerance(1e-6);
    mystepper.SetMaxiters(10);
    mystepper.SetMaxItersSuccess(3);
    mystepper.SetRequiredSuccessfulSteps(5);
    mystepper.SetStepIncreaseFactor(2);
    mystepper.SetStepDecreaseFactor(0.5);
    mystepper.SetVerbose(true);

    // Execute the time integration
    for (int i = 0; i < num_steps; i++) {
        mystepper.Advance(timestep);
        printf("    %7.4f | %12.8f  %12.8f  %12.8f  %12.8f | %12.8f | %2d\n",
               mystepper.GetTime(),                         // time
               mystepper.get_X()(0), mystepper.get_X()(1),  // position
               mystepper.get_V()(0), mystepper.get_V()(1),  // velocity
               mystepper.get_L()(0),                        // Lagrange multiplier
               mystepper.GetNumIterations()                 // number of HHT iterations
               );
    }
}

// ==========================================================================================================

int main(int argc, char* argv[]) {
    //Oscillator();
    Pendulum();
    return 0;
}
