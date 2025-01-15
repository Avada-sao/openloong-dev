#include <mujoco/mujoco.h>
#include <GLFW/glfw3.h>
#include <cstdio>
#include <iostream>
#include "useful_math.h"
#include "GLFW_callbacks.h"
#include "MJ_interface.h"
#include "PVT_ctrl.h"
#include "pino_kin_dyn.h"
#include "data_logger.h"
#include "wbc_priority.h"
#include "gait_scheduler.h"
#include "foot_placement.h"
#include "joystick_interpreter.h"

// MuJoCo load and compile model
char error[1000] = "Could not load binary model";
mjModel* mj_model = mj_loadXML("../models/scene.xml", 0, error, 1000);
mjData* mj_data = mj_makeData(mj_model);

//************************
int main(int argc, const char** argv){
    // ini classes
    UIctr uiController(mj_model,mj_data);   // UI control for Mujoco
    MJ_Interface mj_interface(mj_model, mj_data); // data interface for Mujoco
    Pin_KinDyn kinDynSolver("../models/centaur.urdf"); // kinematics and dynamics solver
    DataBus RobotState(kinDynSolver.model_nv); // data bus
    WBC_priority WBC_solv(kinDynSolver.model_nv, 18, 26, 0.5, mj_model->opt.timestep); // WBC solver
    GaitScheduler gaitScheduler(0.25, mj_model->opt.timestep); // gait scheduler
    PVT_Ctr pvtCtr(mj_model->opt.timestep,"../common/joint_ctrl_config.json");// PVT joint control
    FootPlacement footPlacement; // foot-placement planner
    JoyStickInterpreter jsInterp(mj_model->opt.timestep); // desired baselink velocity generator
    DataLogger logger("../record/datalog.log"); // data logger

    // variables ini
    double stand_legLength = 0.43; // desired baselink height
    double foot_height = 0.07; // distance between the foot ankel joint and the bottom
    double xv_des = 0.5;  // desired velocity in x direction

    RobotState.width_hips = 0.332;
    footPlacement.kp_vx = 0.03;
    footPlacement.kp_vy = 0.05;
    footPlacement.kp_wz = 0.0;
    footPlacement.stepHeight = 0.15;
    footPlacement.legLength=stand_legLength;
    //mju_copy(mj_data->qpos, mj_model->key_qpos, mj_model->nq*1); // set ini pos in Mujoco
    int model_nv=kinDynSolver.model_nv;

    // ini position and posture for foot-end and hand
    std::vector<double> motors_pos_des(model_nv-6,0);
    std::vector<double> motors_pos_cur(model_nv-6,0);
    std::vector<double> motors_vel_des(model_nv-6,0);
    std::vector<double> motors_vel_cur(model_nv-6,0);
    std::vector<double> motors_tau_des(model_nv-6,0);
    std::vector<double> motors_tau_cur(model_nv-6,0);
    
    // register variable name for data logger
    logger.addIterm("simTime", 1);
    logger.addIterm("motors_pos_cur",model_nv-6);
    logger.addIterm("motors_vel_cur",model_nv-6);
    logger.addIterm("rpy",3);
    logger.addIterm("fL",3);
    logger.addIterm("fR",3);
    logger.addIterm("basePos",3);
    logger.addIterm("baseLinVel",3);
    logger.addIterm("baseAcc",3);
    logger.finishItermAdding();

    /// ----------------- sim Loop ---------------
    double simEndTime=60;
    mjtNum simstart = mj_data->time;
    double simTime = mj_data->time;
    double openLoopCtrTime=3;
    double startSteppingTime=7;
    double startWalkingTime=10;

    // init UI: GLFW
    uiController.iniGLFW();
    uiController.enableTracking(); // enable viewpoint tracking of the body 1 of the robot
    uiController.createWindow("Demo",false);
    UIctr::ButtonState buttonState;

    while(!glfwWindowShouldClose(uiController.window)){
        // advance interactive simulation for 1/60 sec
        //  Assuming MuJoCo can simulate faster than real-time, which it usually can,
        //  this loop will finish on time for the next frame to be rendered at 60 fps.
        //  Otherwise add a cpu timer and exit this loop when it is time to render.
        simstart=mj_data->time;
        while(mj_data->time - simstart < 1.0/60.0 && uiController.runSim){
            mj_step(mj_model, mj_data);

            simTime=mj_data->time;
            printf("-------------%.3f s------------\n",simTime);
            mj_interface.updateSensorValues();
            mj_interface.dataBusWrite(RobotState);

            // input from joystick
            // space: start and stop stepping (after 3s)
            // w: forward walking
            // s: stop forward walking
            // a: turning left
            // d: turning right
            buttonState=uiController.getButtonState();
            std::cout<<"motionState: "<<RobotState.motionState<<std::endl;
            if (simTime > openLoopCtrTime){
                // if (buttonState.key_space && RobotState.motionState==DataBus::Stand) {
                //     jsInterp.setIniPos(RobotState.q(0), RobotState.q(1), RobotState.base_rpy(2));
                    RobotState.motionState = DataBus::Walk;
                //     std::cout<<buttonState.key_space<<std::endl;
                // }else if (buttonState.key_space && RobotState.motionState==DataBus::Walk && fabs(jsInterp.vxLGen.y)<0.01) {
                //     // RobotState.motionState = DataBus::Walk2Stand;
                //     jsInterp.setIniPos(RobotState.q(0), RobotState.q(1), RobotState.base_rpy(2));
                //     std::cout<<buttonState.key_space<<std::endl;
                // }

                // if (buttonState.key_a && RobotState.motionState!=DataBus::Stand) {
                //     if (jsInterp.wzLGen.yDes<0)
                //         jsInterp.setWzDesLPara(0, 0.5);
                //     else
                        // jsInterp.setWzDesLPara(-0.1, 1.0);
                // }
                // if (buttonState.key_d && RobotState.motionState!=DataBus::Stand) {
                //     if (jsInterp.wzLGen.yDes>0)
                //         jsInterp.setWzDesLPara(0, 0.5);
                //     else
                        jsInterp.setWzDesLPara(-0.35, 2.0);
                // }

                // if (buttonState.key_w && RobotState.motionState!=DataBus::Stand)
                    // jsInterp.setVyDesLPara(0.3, 1.0);

                // if (buttonState.key_s && RobotState.motionState!=DataBus::Stand)
                    // jsInterp.setVxDesLPara(0, 0.5);

                // if (buttonState.key_h)
                //     jsInterp.setIniPos(RobotState.q(0), RobotState.q(1), RobotState.base_rpy(2));
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

            Eigen::Vector3d xyz;

            // switch between walk and stand
            if (RobotState.motionState==DataBus::Walk || RobotState.motionState==DataBus::Walk2Stand) {
                jsInterp.step();
                RobotState.js_pos_des(2) = stand_legLength; // pos z is not assigned in jyInterp
                jsInterp.dataBusWrite(RobotState); // only pos x, pos y, theta z, vel x, vel y , omega z are rewrote.

//                if (simTime <startSteppingTime+0.002)
//                    RobotState.motionState=DataBus::Walk;
                // gait scheduler
                gaitScheduler.dataBusRead(RobotState);
                gaitScheduler.step();
                gaitScheduler.dataBusWrite(RobotState);

                footPlacement.dataBusRead(RobotState);
                footPlacement.getSwingPos();
                footPlacement.dataBusWrite(RobotState);

                xyz(0) = RobotState.js_pos_des(0);
                xyz(1) = RobotState.js_pos_des(1);
                xyz(2) = stand_legLength;
            }
            else{
                xyz = RobotState.q.block(0,0,3,1);
            }

            // if (simTime <= openLoopCtrTime || RobotState.motionState==DataBus::Walk2Stand) {
            //     WBC_solv.setQini(initJointPos, RobotState.q);
            //     WBC_solv.fe_fl_pos_des_W=RobotState.fe_fl_pos_W;
            //     WBC_solv.fe_fr_pos_des_W=RobotState.fe_fr_pos_W;
            //     WBC_solv.fe_rl_pos_des_W=RobotState.fe_rl_pos_W;
            //     WBC_solv.fe_rr_pos_des_W=RobotState.fe_rr_pos_W;
            //     WBC_solv.fe_fl_rot_des_W=RobotState.fe_fl_rot_W;
            //     WBC_solv.fe_fr_rot_des_W=RobotState.fe_fr_rot_W;
            //     WBC_solv.fe_rl_rot_des_W=RobotState.fe_rl_rot_W;
            //     WBC_solv.fe_rr_rot_des_W=RobotState.fe_rr_rot_W;
            //     WBC_solv.pCoMDes= RobotState.pCoM_W;
            //     WBC_solv.pCoMDes(0)=(RobotState.fe_fl_pos_W(0)+RobotState.fe_fr_pos_W(0))*0.5;
            //     WBC_solv.pCoMDes(1)=(RobotState.fe_fl_pos_W(1)+RobotState.fe_fr_pos_W(1))*0.5;
            // }

            // if (RobotState.motionState==DataBus::Stand) {
            //     WBC_solv.pCoMDes(0) = (RobotState.fe_fl_pos_W(0) + RobotState.fe_fr_pos_W(0)) * 0.5;
            //     WBC_solv.pCoMDes(1) = (RobotState.fe_fl_pos_W(1) + RobotState.fe_fr_pos_W(1)) * 0.5;
            // }

//            std::cout<<"pCoM_W"<<std::endl<<RobotState.pCoM_W.transpose()<<std::endl<<"pCoM_Des"<<std::endl<<WBC_solv.pCoMDes.transpose()<<std::endl;

            // ------------- WBC ------------
            // WBC input
            RobotState.Fr_ff = Eigen::VectorXd::Zero(12);
            RobotState.des_ddq = Eigen::VectorXd::Zero(mj_model->nv);
            RobotState.des_dq = Eigen::VectorXd::Zero(mj_model->nv);
            RobotState.des_delta_q = Eigen::VectorXd::Zero(mj_model->nv);
            RobotState.base_rpy_des << 0, 0, jsInterp.thetaZ;
            RobotState.base_pos_des = xyz;

            RobotState.Fr_ff<<0,0,70,  0,0,70,  0,0,70,  0,0,70;

            // adjust des_delata_q, des_dq and des_ddq to achieve forward walking
            if (RobotState.motionState==DataBus::Walk) {
                RobotState.des_delta_q.block<2, 1>(0, 0) << jsInterp.vx_W * mj_model->opt.timestep, jsInterp.vy_W * mj_model->opt.timestep;
                RobotState.des_delta_q(5) = jsInterp.wz_L * mj_model->opt.timestep;
                RobotState.des_dq.block<2, 1>(0, 0) << jsInterp.vx_W, jsInterp.vy_W;
                RobotState.des_dq(5) = jsInterp.wz_L;

                double k = 5; //5
                RobotState.des_ddq.block<2, 1>(0, 0) << k * (jsInterp.vx_W - RobotState.dq(0)), k * (jsInterp.vy_W -
                                                                                                     RobotState.dq(1));
                RobotState.des_ddq(5) = k * (jsInterp.wz_L - RobotState.dq(5));
            }
            printf("js_vx=%.3f js_vy=%.3f wz_L=%.3f px_w=%.3f py_w=%.3f thetaZ=%.3f\n", jsInterp.vx_W,jsInterp.vy_W,jsInterp.wz_L, jsInterp.px_W, jsInterp.py_W, jsInterp.thetaZ);

            // WBC Calculation
            WBC_solv.dataBusRead(RobotState);
            WBC_solv.computeDdq(kinDynSolver);
            WBC_solv.computeTau();
            WBC_solv.dataBusWrite(RobotState);

            Eigen::VectorXd initJointPos=Eigen::VectorXd::Zero(12);
            initJointPos << 0.0, 0.75, -1.5, 0.0, 0.75, -1.5, 0.0, 0.75, -1.5, 0.0, 0.75, -1.5;

            // get the final joint command
            if (simTime<=openLoopCtrTime || RobotState.motionState==DataBus::Stand){
                RobotState.motors_pos_des= eigen2std(initJointPos);
                RobotState.motors_vel_des=motors_vel_des;
                RobotState.motors_tor_des=motors_tau_des;
                std::cout<<"huohuohuohuohuohuouhuo"<<std::endl;
            }else{
                Eigen::VectorXd pos_des=kinDynSolver.integrateDIY(RobotState.q, RobotState.wbc_delta_q_final);
                RobotState.motors_pos_des= eigen2std(pos_des.block(7,0, model_nv-6,1));
                RobotState.motors_vel_des = eigen2std(RobotState.wbc_dq_final);
                RobotState.motors_tor_des = eigen2std(RobotState.wbc_tauJointRes);
            }

            pvtCtr.dataBusRead(RobotState);
            if (simTime<=openLoopCtrTime){
                pvtCtr.calMotorsPVT(100.0/1000.0/180.0*3.1415);
            }else{
                if (RobotState.motionState == DataBus::Walk2Stand || RobotState.motionState == DataBus::Walk){
                    pvtCtr.setJointPD(400,20,"FL_hip_joint");
                    pvtCtr.setJointPD(400,20,"RR_hip_joint");
                    pvtCtr.setJointPD(300,10,"FR_hip_joint");
                    pvtCtr.setJointPD(400,20,"RL_hip_joint");
                    pvtCtr.setJointPD(400,20,"FL_thigh_joint");
                    pvtCtr.setJointPD(300,10,"FR_thigh_joint");
                    pvtCtr.setJointPD(400,20,"RR_thigh_joint");
                    pvtCtr.setJointPD(400,20,"RL_thigh_joint");
                    pvtCtr.setJointPD(400,10,"FL_calf_joint");
                    pvtCtr.setJointPD(300,5,"FR_calf_joint");
                    pvtCtr.setJointPD(400,10,"RR_calf_joint");
                    pvtCtr.setJointPD(400,10,"RL_calf_joint");
                }else{
                    pvtCtr.setJointPD(200,100,"FL_hip_joint");
                    pvtCtr.setJointPD(200,100,"RR_hip_joint");
                    pvtCtr.setJointPD(200,100,"FR_hip_joint");
                    pvtCtr.setJointPD(200,100,"RL_hip_joint");
                    pvtCtr.setJointPD(400,20,"FL_thigh_joint");
                    pvtCtr.setJointPD(400,20,"FR_thigh_joint");
                    pvtCtr.setJointPD(400,20,"RR_thigh_joint");
                    pvtCtr.setJointPD(400,20,"RL_thigh_joint");
                    pvtCtr.setJointPD(400,10,"FL_calf_joint");
                    pvtCtr.setJointPD(400,10,"FR_calf_joint");
                    pvtCtr.setJointPD(400,10,"RR_calf_joint");
                    pvtCtr.setJointPD(400,10,"RL_calf_joint");
                }
                pvtCtr.calMotorsPVT();
            }
            pvtCtr.dataBusWrite(RobotState);
            mj_interface.setMotorsTorque(RobotState.motors_tor_out);

            // data record` 1111

            logger.startNewLine();
            logger.recItermData("simTime", simTime);
            logger.recItermData("motors_pos_cur",RobotState.motors_pos_cur);
            logger.recItermData("motors_vel_cur",RobotState.motors_vel_cur);
            logger.recItermData("rpy",RobotState.rpy);
            logger.recItermData("fL",RobotState.fL);
            logger.recItermData("fR",RobotState.fR);
            logger.recItermData("basePos",RobotState.basePos);
            logger.recItermData("baseLinVel",RobotState.baseLinVel);
            logger.recItermData("baseAcc",RobotState.baseAcc);
            logger.finishLine();

            printf("rpyVal=[%.5f, %.5f, %.5f]\n", RobotState.rpy[0], RobotState.rpy[1], RobotState.rpy[2]);
            printf("gps=[%.5f, %.5f, %.5f]\n", RobotState.basePos[0], RobotState.basePos[1], RobotState.basePos[2]);
            printf("vel=[%.5f, %.5f, %.5f]\n", RobotState.baseLinVel[0], RobotState.baseLinVel[1], RobotState.baseLinVel[2]);
        }

        if(mj_data->time>=simEndTime){
            break;
        }
        uiController.updateScene();
    }

//    // free visualization storage
    uiController.Close();

    return 0;
}