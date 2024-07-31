#include <mujoco/mujoco.h>
#include <GLFW/glfw3.h>
#include "GLFW_callbacks.h"
#include "MJ_interface.h"
#include "PVT_ctrl.h"
#include "data_logger.h"
#include "data_bus.h"
#include "pino_kin_dyn.h"
#include "useful_math.h"
#include "wbc_priority.h"
#include "mpc.h"
#include "gait_scheduler.h"
#include "foot_placement.h"
#include "joystick_interpreter.h"
#include <string>
#include <iostream>

const   double  dt = 0.001;
const   double  dt_200Hz = 0.005;
// MuJoCo load and compile model
char error[1000] = "Could not load binary model";
mjModel* mj_model = mj_loadXML("../models/scene.xml", 0, error, 1000);
mjData* mj_data = mj_makeData(mj_model);

int main(int argc, char **argv) {
    // initialize classes
    UIctr uiController(mj_model,mj_data);   // UI control for Mujoco
    MJ_Interface mj_interface(mj_model, mj_data); // data interface for Mujoco
    Pin_KinDyn kinDynSolver("../models/AzureLoong.urdf"); // kinematics and dynamics solver
    DataBus RobotState(kinDynSolver.model_nv); // data bus
    WBC_priority WBC_solv(kinDynSolver.model_nv, 18, 22, 0.7, mj_model->opt.timestep); // WBC solver
    MPC MPC_solv(dt_200Hz);  // mpc controller
    GaitScheduler gaitScheduler(0.3, mj_model->opt.timestep); // gait scheduler
    PVT_Ctr pvtCtr(mj_model->opt.timestep,"../common/joint_ctrl_config.json");// PVT joint control
    FootPlacement footPlacement; // foot-placement planner
    JoyStickInterpreter jsInterp(mj_model->opt.timestep); // desired baselink velocity generator
    DataLogger logger("../record/datalog.log"); // data logger

    // initialize UI: GLFW
    uiController.iniGLFW();
    uiController.enableTracking(); // enable viewpoint tracking of the body 1 of the robot
    uiController.createWindow("Demo",false); // NOTE: if the saveVideo is set to true, the raw recorded file could be 2.5 GB for 15 seconds!
	UIctr::ButtonState buttonState;

    // initialize variables
    double stand_legLength = 1.01;//-0.95; // desired baselink height
    double foot_height =0.07; // distance between the foot ankel joint and the bottom
    double xv_des = 1.2;  // desired velocity in x direction

    const int robot_nq= kinDynSolver.model_nv + 1;
    const int robot_nv= robot_nq + 1;

    RobotState.width_hips = 0.229;
    footPlacement.kp_vx = 0.03;
    footPlacement.kp_vy = 0.03;
    footPlacement.kp_wz = 0.03;
    footPlacement.stepHeight = 0.15;
    footPlacement.legLength=stand_legLength;

    mju_copy(mj_data->qpos, mj_model->key_qpos, mj_model->nq*1); // set ini pos in Mujoco

    std::vector<double> motors_pos_des(robot_nv - 6, 0);
    std::vector<double> motors_pos_cur(robot_nv - 6, 0);
    std::vector<double> motors_vel_des(robot_nv - 6, 0);
    std::vector<double> motors_vel_cur(robot_nv - 6, 0);
    std::vector<double> motors_tau_des(robot_nv - 6, 0);
    std::vector<double> motors_tau_cur(robot_nv - 6, 0);

    // ini position and posture for foot-end and hand
    Eigen::Vector3d fe_l_pos_L_des={-0.018, 0.113, -stand_legLength};
    Eigen::Vector3d fe_r_pos_L_des={-0.018, -0.116, -stand_legLength};
    Eigen::Vector3d fe_l_eul_L_des={-0.000, -0.008, -0.000};
    Eigen::Vector3d fe_r_eul_L_des={0.000, -0.008, 0.000};
    Eigen::Matrix3d fe_l_rot_des= eul2Rot(fe_l_eul_L_des(0),fe_l_eul_L_des(1),fe_l_eul_L_des(2));
    Eigen::Matrix3d fe_r_rot_des= eul2Rot(fe_r_eul_L_des(0),fe_r_eul_L_des(1),fe_r_eul_L_des(2));

    Eigen::Vector3d hd_l_pos_L_des={-0.02, 0.32, -0.159};
    Eigen::Vector3d hd_r_pos_L_des={-0.02, -0.32, -0.159};
    Eigen::Vector3d hd_l_eul_L_des={-1.7581, 0.2129, 2.9581};
    Eigen::Vector3d hd_r_eul_L_des={1.7581, 0.2129, -2.9581};
    Eigen::Matrix3d hd_l_rot_des= eul2Rot(hd_l_eul_L_des(0),hd_l_eul_L_des(1),hd_l_eul_L_des(2));
    Eigen::Matrix3d hd_r_rot_des= eul2Rot(hd_r_eul_L_des(0),hd_r_eul_L_des(1),hd_r_eul_L_des(2));

    auto resLeg=kinDynSolver.computeInK_Leg(fe_l_rot_des,fe_l_pos_L_des,fe_r_rot_des,fe_r_pos_L_des);
    auto resHand=kinDynSolver.computeInK_Hand(hd_l_rot_des,hd_l_pos_L_des,hd_r_rot_des,hd_r_pos_L_des);
    Eigen::VectorXd qIniDes=Eigen::VectorXd::Zero(robot_nq, 1);
    qIniDes.block(7, 0, robot_nq - 7, 1)= resLeg.jointPosRes + resHand.jointPosRes;

    // register variable name for data logger
    logger.addIterm("simTime",1);
    logger.addIterm("motor_pos_des",30);
    logger.addIterm("motor_pos_cur",30);
    logger.addIterm("motor_vel_des",30);
    logger.addIterm("motor_vel_cur",30);
    logger.addIterm("motor_tor_des",30);
    logger.addIterm("rpyVal",3);
    logger.addIterm("base_omega_W",3);
    logger.addIterm("gpsVal",3);
    logger.addIterm("base_vel",3);
	logger.addIterm("dX_cal",12);
	logger.addIterm("Ufe",12);
    logger.finishItermAdding();

    //// -------------------------- main loop --------------------------------

    int  MPC_count = 0; // count for controlling the mpc running period

    double openLoopCtrTime=3;
    double startSteppingTime=7;
    double startWalkingTime=10;
    double simEndTime=200;

    mjtNum simstart = mj_data->time;
    double simTime = mj_data->time;

    while (!glfwWindowShouldClose(uiController.window)) {
        simstart = mj_data->time;
        while (mj_data->time - simstart < 1.0 / 60.0 && uiController.runSim) {
            mj_step(mj_model, mj_data);
            simTime=mj_data->time;
            // Read the sensors:
            mj_interface.updateSensorValues();
            mj_interface.dataBusWrite(RobotState);

			// input from joystick
			// space: start and stop stepping (after 3s)
			// w: forward walking
			// s: stop forward walking
			// a: turning left
			// d: turning right
			buttonState=uiController.getButtonState();
			if (simTime > openLoopCtrTime){
				if (buttonState.key_space && RobotState.motionState==DataBus::Stand) {
					jsInterp.setIniPos(RobotState.q(0), RobotState.q(1), RobotState.base_rpy(2));
					RobotState.motionState = DataBus::Walk;
				}else if (buttonState.key_space && RobotState.motionState==DataBus::Walk && fabs(jsInterp.vxLGen.y)<0.01) {
					RobotState.motionState = DataBus::Walk2Stand;
					jsInterp.setIniPos(RobotState.q(0), RobotState.q(1), RobotState.base_rpy(2));
				}

				if (buttonState.key_a && RobotState.motionState!=DataBus::Stand) {
					if (jsInterp.wzLGen.yDes<0)
						jsInterp.setWzDesLPara(0, 0.5);
					else
						jsInterp.setWzDesLPara(0.35, 1.0);
				}
				if (buttonState.key_d && RobotState.motionState!=DataBus::Stand) {
					if (jsInterp.wzLGen.yDes>0)
						jsInterp.setWzDesLPara(0, 0.5);
					else
						jsInterp.setWzDesLPara(-0.35, 1.0);
				}

				if (buttonState.key_w && RobotState.motionState!=DataBus::Stand)
					jsInterp.setVxDesLPara(xv_des, 2.0);

				if (buttonState.key_s && RobotState.motionState!=DataBus::Stand)
					jsInterp.setVxDesLPara(0, 0.5);

				if (buttonState.key_h)
					jsInterp.setIniPos(RobotState.q(0), RobotState.q(1), RobotState.base_rpy(2));
			}

            // update kinematics and dynamics info
            kinDynSolver.dataBusRead(RobotState);
            kinDynSolver.computeJ_dJ();
            kinDynSolver.computeDyn();
            kinDynSolver.dataBusWrite(RobotState);

			if (simTime>=openLoopCtrTime && simTime<openLoopCtrTime+0.002) {
				RobotState.motionState = DataBus::Stand;
			}

			if (RobotState.motionState==DataBus::Walk2Stand || simTime<= openLoopCtrTime)
				jsInterp.setIniPos(RobotState.q(0), RobotState.q(1), RobotState.base_rpy(2));


			// switch between walk and stand
			if (RobotState.motionState==DataBus::Walk || RobotState.motionState==DataBus::Walk2Stand) {
				jsInterp.step();
				RobotState.js_pos_des(2) = stand_legLength + foot_height; // pos z is not assigned in jyInterp
				jsInterp.dataBusWrite(RobotState); // only pos x, pos y, theta z, vel x, vel y , omega z are rewrote.

                MPC_solv.enable();

                // gait scheduler
                gaitScheduler.dataBusRead(RobotState);
                gaitScheduler.step();
                gaitScheduler.dataBusWrite(RobotState);

                footPlacement.dataBusRead(RobotState);
                footPlacement.getSwingPos();
                footPlacement.dataBusWrite(RobotState);
			}

			if (simTime <= openLoopCtrTime || RobotState.motionState==DataBus::Walk2Stand) {
                WBC_solv.setQini(qIniDes, RobotState.q);
                WBC_solv.fe_l_pos_des_W=RobotState.fe_l_pos_W;
                WBC_solv.fe_r_pos_des_W=RobotState.fe_r_pos_W;
                WBC_solv.fe_l_rot_des_W=RobotState.fe_l_rot_W;
                WBC_solv.fe_r_rot_des_W=RobotState.fe_r_rot_W;
                WBC_solv.pCoMDes= RobotState.pCoM_W;
            }

            // ------------- MPC ------------
			MPC_count = MPC_count + 1;
            if (MPC_count > (dt_200Hz / dt - 1)) { //MPC_count = 1, 2, 3, 4, 5(5 run MPC)
                MPC_solv.dataBusRead(RobotState);
                MPC_solv.cal();
                MPC_solv.dataBusWrite(RobotState);
                MPC_count = 0;
            }
            // ------------- WBC ------------
            // WBC Calculation
            WBC_solv.dataBusRead(RobotState);
            WBC_solv.computeDdq(kinDynSolver);
            WBC_solv.computeTau();
            WBC_solv.dataBusWrite(RobotState);

            // get the final joint command
            if (simTime <= openLoopCtrTime) {
                RobotState.motors_pos_des = eigen2std(resLeg.jointPosRes + resHand.jointPosRes);
                RobotState.motors_vel_des = motors_vel_des;
                RobotState.motors_tor_des = motors_tau_des;
            } else {

                Eigen::Matrix<double, 1, nx>  L_diag;
                Eigen::Matrix<double, 1, nu>  K_diag;
                L_diag <<
                       1.0, 1.0, 1.0,//eul
                        1.0, 200.0,  1.0,//pCoM
                        1e-7, 1e-7, 1e-7,//w
                        100.0, 100.0, 1.0;//vCoM
                K_diag <<
                       1.0, 1.0, 1.0,//fl
                        1.0, 1.0, 1.0,
                        1.0, 1.0, 1.0,//fr
                        1.0, 1.0, 1.0,1.0;
                MPC_solv.set_weight(1e-6, L_diag, K_diag);

                Eigen::VectorXd pos_des = kinDynSolver.integrateDIY(RobotState.q, RobotState.wbc_delta_q_final);
                RobotState.motors_pos_des = eigen2std(pos_des.block(7, 0, robot_nv - 6, 1));
                RobotState.motors_vel_des = eigen2std(RobotState.wbc_dq_final);
                RobotState.motors_tor_des = eigen2std(RobotState.wbc_tauJointRes);
            }

            // joint PVT controller
            pvtCtr.dataBusRead(RobotState);
            if (simTime <= openLoopCtrTime) {
                pvtCtr.calMotorsPVT(110.0 / 1000.0 / 180.0 * 3.1415);
            } else {
                pvtCtr.setJointPD(100,10,"J_ankle_l_pitch");
                pvtCtr.setJointPD(100,10,"J_ankle_l_roll");
                pvtCtr.setJointPD(100,10,"J_ankle_r_pitch");
                pvtCtr.setJointPD(100,10,"J_ankle_r_roll");
                pvtCtr.setJointPD(1000,100,"J_knee_l_pitch");
                pvtCtr.setJointPD(1000,100,"J_knee_r_pitch");
                pvtCtr.calMotorsPVT();
            }
            pvtCtr.dataBusWrite(RobotState);

            // give the joint torque command to Webots
            mj_interface.setMotorsTorque(RobotState.motors_tor_out);

            // data save
            logger.startNewLine();
            logger.recItermData("simTime", simTime);
            logger.recItermData("motor_pos_des", RobotState.motors_pos_des);
            logger.recItermData("motor_pos_cur", RobotState.motors_pos_cur);
            logger.recItermData("motor_vel_des", RobotState.motors_vel_des);
            logger.recItermData("motor_vel_cur", RobotState.motors_vel_cur);
            logger.recItermData("motor_tor_des", RobotState.motors_tor_des);
            logger.recItermData("rpyVal", RobotState.rpy);
            logger.recItermData("base_omega_W", RobotState.base_omega_W);
            logger.recItermData("gpsVal", RobotState.basePos);
            logger.recItermData("base_vel", RobotState.dq.block<3, 1>(0, 0));
			logger.recItermData("dX_cal",RobotState.dX_cal);
			logger.recItermData("Ufe",RobotState.Fr_ff);
            logger.finishLine();
        }

        if (mj_data->time>=simEndTime)
            break;

        uiController.updateScene();
    };
    // free visualization storage
    uiController.Close();

    return 0;
}
