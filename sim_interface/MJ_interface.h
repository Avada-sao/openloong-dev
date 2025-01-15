/*
This is part of OpenLoong Dynamics Control, an open project for the control of biped robot,
Copyright (C) 2024 Humanoid Robot (Shanghai) Co., Ltd, under Apache 2.0.
Feel free to use in any purpose, and cite OpenLoong-Dynamics-Control in any style, to contribute to the advancement of the community.
 <https://atomgit.com/openloong/openloong-dyn-control.git>
 <web@openloong.org.cn>
*/
#pragma once

#include <mujoco/mujoco.h>
#include "data_bus.h"
#include <string>
#include <vector>

class MJ_Interface {
public:
    int jointNum{0};
    std::vector<double> motor_pos;
    std::vector<double> motor_pos_Old;
    std::vector<double> motor_vel;
    double rpy[3]{0}; // roll,pitch and yaw of baselink
    double baseQuat[4]{0}; // in quat, mujoco order is [w,x,y,z], here we rearrange to [x,y,z,w]
    double f3d[3][2]{0}; // 3D foot-end contact force, L for 1st col, R for 2nd col
    double basePos[3]{0}; // position of baselink, in world frame
    double baseAcc[3]{0};  // acceleration of baselink, in body frame
    double baseAngVel[3]{0}; // angular velocity of baselink, in body frame
    double baseLinVel[3]{0}; // linear velocity of baselink, in body frame
    const std::vector<std::string> JointName={"FL_hip_joint", "FL_thigh_joint", "FL_calf_joint",
                                              "FR_hip_joint", "FR_thigh_joint", "FR_calf_joint",
                                              "RL_hip_joint", "RL_thigh_joint", "RL_calf_joint",
                                              "RR_hip_joint", "RR_thigh_joint", "RR_calf_joint"}; // joint name in XML file, the corresponds motors name should be M_*, ref to line 29 of MJ_Interface.cpp
    const std::string baseName="base";
    const std::string orientationSensorName="baselink-quat"; // in quat, mujoco order is [w,x,y,z], here we rearrange to [x,y,z,w]
    const std::string velSensorName="baselink-velocity";
    const std::string gyroSensorName="baselink-gyro";
    const std::string accSensorName="baselink-baseAcc";

    MJ_Interface(mjModel *mj_modelIn, mjData  *mj_dataIn);
    void updateSensorValues();
    void setMotorsTorque(std::vector<double> &tauIn);
    void dataBusWrite(DataBus &busIn);

private:
    mjModel *mj_model;
    mjData  *mj_data;
    std::vector<int> jntId_qpos, jntId_qvel, jntId_dctl;

    int orientataionSensorId;
    int velSensorId;
    int gyroSensorId;
    int accSensorId;
    int baseBodyId;

    double timeStep{0.001}; // second
    bool isIni{false};
};



