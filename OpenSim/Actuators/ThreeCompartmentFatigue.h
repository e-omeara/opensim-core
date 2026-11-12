#ifndef OPENSIM_THREE_COMPARTMENT_FATIGUE_H_
#define OPENSIM_THREE_COMPARTMENT_FATIGUE_H_
/* -------------------------------------------------------------------------- *
 *                     OpenSim:  ThreeCompartmentFatigue.h                    *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include "osimActuatorsDLL.h"
#include <OpenSim/Simulation/Model/ModelComponent.h>
#include <OpenSim/Simulation/Model/Muscle.h>

#ifdef SWIG
    #ifdef OSIMACTUATORS_API
        #undef OSIMACTUATORS_API
        #define OSIMACTUATORS_API
    #endif
#endif

namespace OpenSim {

//==============================================================================
//                       THREE COMPARTMENT FATIGUE
//==============================================================================
/** 
 * TODO: [add description here]
*/
class OSIMACTUATORS_API ThreeCompartmentFatigue : public ModelComponent {
OpenSim_DECLARE_CONCRETE_OBJECT(ThreeCompartmentFatigue, ModelComponent);

public:
//==============================================================================
// PROPERTIES
//==============================================================================
    OpenSim_DECLARE_PROPERTY(fatigue_rate, double,
        "Fatigue rate F (1/s): rate at which active units become fatigued.");
    OpenSim_DECLARE_PROPERTY(recovery_rate, double,
        "Recovery rate R (1/s): rate at which fatigued units recover.");
    OpenSim_DECLARE_PROPERTY(rest_recovery_multiplier, double,
        "Multiplier r applied to the recovery rate while the muscle rests "
        "(3CC-r). Use 1 for the original 3CC model.");
    OpenSim_DECLARE_PROPERTY(development_factor, double,
        "L_D (1/s): rate at which resting units are recruited when the target "
        "load is above the active fraction.");
    OpenSim_DECLARE_PROPERTY(relaxation_factor, double,
        "L_R (1/s): rate at which active units are released when the target "
        "load is below the active fraction.");
    OpenSim_DECLARE_PROPERTY(rest_target_load, double,
        "The muscle is resting (recovery multiplied by r) when the target "
        "load is at or below this value. Default 0, which never occurs while "
        "a muscle has non-zero activation.");
    OpenSim_DECLARE_PROPERTY(initial_active_fraction, double,
        "Initial value of M_A.");
    OpenSim_DECLARE_PROPERTY(initial_fatigued_fraction, double,
        "Initial value of M_F.");

//==============================================================================
// SOCKETS
//==============================================================================
    OpenSim_DECLARE_SOCKET(muscle, Muscle,
        "The muscle whose activation drives this fatigue model.");

//==============================================================================
// OUTPUTS
//==============================================================================
    OpenSim_DECLARE_OUTPUT(force_capacity, double, getForceCapacity,
        SimTK::Stage::Time);
    OpenSim_DECLARE_OUTPUT(fatigued_fraction, double, getFatiguedFraction,
        SimTK::Stage::Time);
    OpenSim_DECLARE_OUTPUT(active_fraction, double, getActiveFraction,
        SimTK::Stage::Time);
    OpenSim_DECLARE_OUTPUT(resting_fraction, double, getRestingFraction,
        SimTK::Stage::Time);

//==============================================================================
// METHODS
//==============================================================================
    ThreeCompartmentFatigue();

    /** Fraction of the maximum isometric force available: 1 - M_F. */
    double getForceCapacity(const SimTK::State& s) const;
    double getFatiguedFraction(const SimTK::State& s) const;
    double getActiveFraction(const SimTK::State& s) const;
    double getRestingFraction(const SimTK::State& s) const;
    void setFatiguedFraction(SimTK::State& s, double value) const;
    void setActiveFraction(SimTK::State& s, double value) const;

    /** The 3CC-r derivatives dM_A/dt and dM_F/dt for the given compartment
    fractions and target load. This is the model's equations without any
    OpenSim state, so they can be tested directly.
    @param[in] activeFraction M_A
    @param[in] fatiguedFraction M_F
    @param[in] targetLoad TL
    @param[out] activeRate dM_A/dt
    @param[out] fatiguedRate dM_F/dt */
    void calcDerivatives(double activeFraction, double fatiguedFraction,
            double targetLoad, double& activeRate, double& fatiguedRate) const;

protected:
    void extendAddToSystem(SimTK::MultibodySystem& system) const override;
    void extendInitStateFromProperties(SimTK::State& s) const override;
    void extendSetPropertiesFromState(const SimTK::State& s) override;
    void computeStateVariableDerivatives(const SimTK::State& s) const override;

private:
    void constructProperties();
};

} // namespace OpenSim

#endif // OPENSIM_THREE_COMPARTMENT_FATIGUE_H_
