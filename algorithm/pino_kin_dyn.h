/*
This is part of OpenLoong Dynamics Control, an open project for the control of biped robot,
Copyright (C) 2024 Humanoid Robot (Shanghai) Co., Ltd, under Apache 2.0.
Feel free to use in any purpose, and cite OpenLoong-Dynamics-Control in any style, to contribute to the advancement of the community.
 <https://atomgit.com/openloong/openloong-dyn-control.git>
 <web@openloong.org.cn>
*/
#pragma once

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/centroidal.hpp"
#include "pinocchio/algorithm/center-of-mass.hpp"
#include "pinocchio/algorithm/aba.hpp"
#include "data_bus.h"
#include <string>
#include "json/json.h"
#include <vector>

class Pin_KinDyn {
public:
    std::vector<bool> motorReachLimit;
    const std::vector<std::string> motorName={"FL_hip_joint", "FL_thigh_joint", "FL_calf_joint",
                                              "FR_hip_joint", "FR_thigh_joint", "FR_calf_joint",
                                              "RL_hip_joint", "RL_thigh_joint", "RL_calf_joint",
                                              "RR_hip_joint", "RR_thigh_joint", "RR_calf_joint"}; // joint name in urdf and jason config files
    Eigen::VectorXd motorMaxTorque;
    Eigen::VectorXd motorMaxPos;
    Eigen::VectorXd motorMinPos;

    Eigen::VectorXd tauJointOld;
    std::string urdf_path;
    pinocchio::Model model_biped;
    pinocchio::Model model_biped_fixed;
    int model_nv;

    pinocchio::FrameIndex fr_foot, fl_foot, rr_foot, rl_foot;
    pinocchio::JointIndex base_joint, fr_hip_joint, fl_hip_joint, rr_hip_joint, rl_hip_joint, fr_thigh_joint, fl_thigh_joint, rr_thigh_joint, rl_thigh_joint;
    pinocchio::JointIndex fr_foot_joint, fl_foot_joint, rr_foot_joint, rl_foot_joint, r_hip_joint_fixed, l_hip_joint_fixed, fr_hip_joint_fixed, fl_hip_joint_fixed, rr_hip_joint_fixed, rl_hip_joint_fixed;
    Eigen::VectorXd q,dq,ddq;
    Eigen::Matrix3d Rcur;
    Eigen::Quaternion<double> quatCur;
    Eigen::Matrix<double,6,-1> J_fr, J_fl, J_rr, J_rl, J_base, J_hip_link;
    Eigen::Matrix<double,6,-1> dJ_fr, dJ_fl, dJ_rr, dJ_rl, dJ_base, dJ_hip_link;
    Eigen::Matrix<double,3,-1> Jcom;
    Eigen::Vector3d fe_fr_pos, fe_fl_pos, fe_rr_pos, fe_rl_pos, base_pos;    // foot-end position in world frame
    Eigen::Vector3d fe_r_pos_body, fe_l_pos_body;  // foot-end position in body frame
    Eigen::Vector3d hd_r_pos, hd_l_pos;  // hand position in world frame
    Eigen::Vector3d hd_r_pos_body, hd_l_pos_body; // hand position in body frame
    Eigen::Vector3d hip_fl_pos, hip_fr_pos, hip_rl_pos, hip_rr_pos, hip_link_pos, thigh_fl_pos, thigh_fr_pos, thigh_rl_pos, thigh_rr_pos;
    Eigen::Vector3d hip_r_pos_body, hip_l_pos_body, hip_fl_pos_body, hip_fr_pos_body, hip_rl_pos_body, hip_rr_pos_body;
    Eigen::Matrix3d hip_link_rot;
    Eigen::Matrix3d fe_r_rot, fe_l_rot, base_rot;
    Eigen::Matrix3d fe_r_rot_body, fe_l_rot_body;
    Eigen::Matrix3d hd_r_rot, hd_l_rot;
    Eigen::Matrix3d hd_r_rot_body, hd_l_rot_body;
    Eigen::MatrixXd legPos;
    bool read;
    
    Eigen::MatrixXd dyn_M, dyn_M_inv, dyn_C, dyn_G, dyn_Ag, dyn_dAg;
    Eigen::VectorXd dyn_Non;
    Eigen::Vector3d CoM_pos;
    Eigen::Matrix3d inertia;
    enum legIdx{
        left,
        right
    };
    struct IkRes{
        int status;
        int itr;
        Eigen::VectorXd err;
        Eigen::VectorXd jointPosRes;
    };

    Pin_KinDyn(std::string urdf_pathIn);
    void dataBusRead(DataBus const &robotState);
    void dataBusWrite(DataBus &robotState);
    void computeJ_dJ();
    Eigen::Vector3d computeFp(Eigen::Vector3d q, int leg);
    void computeDyn();
    IkRes computeInK_Leg(const Eigen::Matrix3d &Rdes_L, const Eigen::Vector3d &Pdes_L, const Eigen::Matrix3d &Rdes_R, const Eigen::Vector3d &Pdes_R);
    IkRes computeInK_Hand(const Eigen::Matrix3d &Rdes_L, const Eigen::Vector3d &Pdes_L, const Eigen::Matrix3d &Rdes_R, const Eigen::Vector3d &Pdes_R);
    Eigen::VectorXd integrateDIY(const Eigen::VectorXd &qI, const Eigen::VectorXd &dqI);
    static Eigen::Quaterniond intQuat(const Eigen::Quaterniond &quat, const Eigen::Matrix<double,3,1> &w);
    void workspaceConstraint(Eigen::VectorXd &qFT, Eigen::VectorXd &tauJointFT);
private:
    pinocchio::Data data_biped, data_biped_fixed;

};
