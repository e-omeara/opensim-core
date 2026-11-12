/* -------------------------------------------------------------------------- *
 *                  OpenSim:  testThreeCompartmentFatigue.cpp                 *
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

#include <OpenSim/Common/Constant.h>
#include <OpenSim/Simulation/osimSimulation.h>
#include <OpenSim/Actuators/osimActuators.h>

#include <catch2/catch_all.hpp>

using namespace OpenSim;

namespace {

// A muscle between ground and a block on a slider, with a fatigue model.
struct FatigueModel {
    Model model;
    Thelen2003Muscle* muscle;
    ThreeCompartmentFatigue* fatigue;

    explicit FatigueModel(bool connectCapacity) {
        auto* block = new Body("block", 1.0, SimTK::Vec3(0),
                SimTK::Inertia(0.1));
        model.addBody(block);
        auto* slider = new SliderJoint("slider", model.getGround(),
                SimTK::Vec3(0), SimTK::Vec3(0), *block, SimTK::Vec3(0),
                SimTK::Vec3(0));
        slider->updCoordinate().setName("x");
        slider->updCoordinate().setDefaultLocked(true);
        model.addJoint(slider);

        muscle = new Thelen2003Muscle("muscle", 100.0, 0.1, 0.2, 0.0);
        muscle->addNewPathPoint("origin", model.getGround(),
                SimTK::Vec3(0, 0.3, 0));
        muscle->addNewPathPoint("insertion", *block, SimTK::Vec3(0));
        model.addForce(muscle);

        fatigue = new ThreeCompartmentFatigue();
        fatigue->setName("fatigue");
        fatigue->connectSocket_muscle(*muscle);
        model.addModelComponent(fatigue);
        if (connectCapacity) {
            muscle->connectInput_force_capacity(
                    fatigue->getOutput("force_capacity"));
        }
        model.finalizeConnections();
    }

    // Drive the muscle with a constant excitation.
    void setConstantExcitation(double excitation) {
        auto* controller = new PrescribedController();
        controller->addActuator(*muscle);
        controller->prescribeControlForActuator("muscle",
                Constant(excitation));
        model.addController(controller);
        model.finalizeConnections();
    }
};

} // namespace

TEST_CASE("ThreeCompartmentFatigue derivative equations") {
    ThreeCompartmentFatigue f;
    double dMA, dMF;

    SECTION("no load, nothing active or fatigued: nothing changes") {
        f.calcDerivatives(0.0, 0.0, 0.0, dMA, dMF);
        CHECK(dMA == 0.0);
        CHECK(dMF == 0.0);
    }

    SECTION("a load above the active fraction recruits resting units") {
        f.calcDerivatives(0.0, 0.0, 0.2, dMA, dMF);
        CHECK(dMA == Catch::Approx(f.get_development_factor() * 0.2));
        CHECK(dMF == 0.0);
    }

    SECTION("recruitment is limited by the resting fraction") {
        // M_R = 0.1, need = 0.8, so C = L_D * M_R.
        f.calcDerivatives(0.2, 0.7, 1.0, dMA, dMF);
        const double expectedC = f.get_development_factor() * 0.1;
        CHECK(dMA == Catch::Approx(expectedC - f.get_fatigue_rate() * 0.2));
    }

    SECTION("a load below the active fraction releases units") {
        f.calcDerivatives(0.5, 0.0, 0.1, dMA, dMF);
        const double expectedC = f.get_relaxation_factor() * (0.1 - 0.5);
        CHECK(dMA == Catch::Approx(expectedC - f.get_fatigue_rate() * 0.5));
        CHECK(dMF == Catch::Approx(f.get_fatigue_rate() * 0.5));
    }

    SECTION("rest recovery is multiplied by r only at rest") {
        f.set_rest_target_load(0.0);
        f.calcDerivatives(0.0, 0.4, 0.0, dMA, dMF);
        CHECK(dMF == Catch::Approx(-f.get_recovery_rate()
                * f.get_rest_recovery_multiplier() * 0.4));
        f.calcDerivatives(0.0, 0.4, 0.01, dMA, dMF);
        CHECK(dMF == Catch::Approx(-f.get_recovery_rate() * 0.4));
    }

    SECTION("the three compartments always sum to 1") {
        // M_R is implicit, so dM_R = -(dM_A + dM_F) by construction; check the
        // explicit form dM_R = -C + R' M_F gives the same rate.
        f.calcDerivatives(0.3, 0.2, 0.6, dMA, dMF);
        const double mr = 1.0 - 0.3 - 0.2;
        const double c = f.get_development_factor() * (0.6 - 0.3);
        const double dMR = -c + f.get_recovery_rate() * 0.2;
        CHECK(-(dMA + dMF) == Catch::Approx(dMR));
        CHECK(mr > 0.6 - 0.3); // the unlimited branch of C applies here
    }
}

TEST_CASE("ThreeCompartmentFatigue state variables and XML round trip") {
    FatigueModel fm(true);
    SimTK::State& s = fm.model.initSystem();

    CHECK(fm.fatigue->getFatiguedFraction(s) == 0.0);
    CHECK(fm.fatigue->getForceCapacity(s) == 1.0);
    CHECK(fm.fatigue->getRestingFraction(s) == 1.0);

    fm.fatigue->setFatiguedFraction(s, 0.25);
    fm.fatigue->setActiveFraction(s, 0.1);
    CHECK(fm.fatigue->getForceCapacity(s) == Catch::Approx(0.75));
    CHECK(fm.fatigue->getRestingFraction(s) == Catch::Approx(0.65));

    fm.fatigue->set_fatigue_rate(0.02);
    fm.model.print("testThreeCompartmentFatigue_model.osim");
    Model copy("testThreeCompartmentFatigue_model.osim");
    const auto& f2 = copy.getComponent<ThreeCompartmentFatigue>("/componentset/fatigue");
    CHECK(f2.get_fatigue_rate() == 0.02);
    copy.initSystem(); // connects the input from its saved path
    CHECK(copy.getComponent<Muscle>("/forceset/muscle")
                  .getInput("force_capacity").isConnected());
}

TEST_CASE("Muscle force_capacity input scales the maximum isometric force") {
    FatigueModel fm(true);
    SimTK::State& s = fm.model.initSystem();

    SECTION("capacity 1 leaves the force unchanged") {
        fm.model.realizeDynamics(s);
        CHECK(fm.muscle->getMaxIsometricForce(s) == Catch::Approx(100.0));
    }

    SECTION("capacity follows the fatigued fraction") {
        fm.fatigue->setFatiguedFraction(s, 0.4);
        fm.model.realizeDynamics(s);
        CHECK(fm.muscle->getMaxIsometricForce(s) == Catch::Approx(60.0));
        // The state-independent accessor still returns the property.
        CHECK(fm.muscle->getMaxIsometricForce() == 100.0);
    }

    SECTION("an unconnected muscle ignores fatigue") {
        FatigueModel unconnected(false);
        SimTK::State& s2 = unconnected.model.initSystem();
        unconnected.fatigue->setFatiguedFraction(s2, 0.4);
        unconnected.model.realizeDynamics(s2);
        CHECK(unconnected.muscle->getMaxIsometricForce(s2)
                == Catch::Approx(100.0));
    }

    SECTION("equilibrium force scales with the capacity") {
        fm.muscle->setActivation(s, 1.0);
        fm.model.equilibrateMuscles(s);
        fm.model.realizeDynamics(s);
        const double fullForce = fm.muscle->getActuation(s);
        REQUIRE(fullForce > 1.0);

        fm.fatigue->setFatiguedFraction(s, 0.5);
        fm.model.equilibrateMuscles(s);
        fm.model.realizeDynamics(s);
        CHECK(fm.muscle->getActuation(s) == Catch::Approx(0.5 * fullForce)
                .epsilon(1e-3));
    }
}

TEST_CASE("ThreeCompartmentFatigue integrates fatigue under constant activation") {
    FatigueModel fm(true);
    // Exaggerate the fatigue rate so a short simulation shows an effect.
    fm.fatigue->set_fatigue_rate(0.1);
    fm.setConstantExcitation(0.5);
    SimTK::State& s = fm.model.initSystem();
    fm.muscle->setActivation(s, 0.5);
    fm.model.equilibrateMuscles(s);

    Manager manager(fm.model);
    manager.initialize(s);
    const SimTK::State& sFinal = manager.integrate(20.0);

    const double mf = fm.fatigue->getFatiguedFraction(sFinal);
    const double ma = fm.fatigue->getActiveFraction(sFinal);
    CHECK(mf > 0.05);
    CHECK(mf < 1.0);
    CHECK(ma > 0.0);
    CHECK(ma + mf <= 1.0);
    CHECK(fm.fatigue->getForceCapacity(sFinal) == Catch::Approx(1.0 - mf));
}
