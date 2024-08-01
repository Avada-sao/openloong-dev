#pragma once

#include "qpOASES.hpp"
#include <algorithm>
#include <Eigen/Dense>
#include "data_bus.h"
#include "useful_math.h"
#include "priority_tasks.h"
#include "pino_kin_dyn.h"
#include <iostream>
#include <iomanip>

class WBC_priority{
public:
    int model_nv;    // size of the system generalized coordinate dq    广义坐标下dp的大小
    Eigen::Vector3d tau_upp_stand_L, tau_low_stand_L; // foot end contact torque limit for stand state, in body frame
    Eigen::Vector3d tau_upp_walk_L, tau_low_walk_L;   // foot end contact torque limit for walk state, in body frame
    double f_z_low{0},f_z_upp{0};                   //Z方向地面反力的最大值与最小值
    DataBus::LegState legStateCur;
    // DataBus::MotionState motionStateCur;
    WBC_priority(int model_nv_In, int QP_nvIn, int QP_ncIn, double miu_In, double dt);      //模型自由度，WBC QP问题的优化变量维数，WBC QP问题的约束维数，模型摩擦系数，控制时域
    double miu{0.5};
    Eigen::MatrixXd dyn_M, dyn_M_inv, dyn_Ag, dyn_dAgdyn_dAg;      //动力学方程参数M，M的逆矩阵，质心动量矩阵
    Eigen::VectorXd dyn_Non;                      // dyn_Non=dyn_C*dq+dyn_G  根据 M*ddq+C*dq+G=tau
    Eigen::MatrixXd Jc, dJc, Jfe, dJfe, Jfe_FL, Jfe_FR, Jfe_HR, Jfe_HL;  
    // 足端接触的雅可比，雅可比的微分，Jfe = [Jfe_L.transpose() Jfe_R.transpose()].transpose() ，左脚右脚的雅可比
    //理论来说Jfe=Jc,但是由于Jfe与足端力矩阵相乘，足端力矩阵大小固定，因此需要固定Jfe大小
    Eigen::MatrixXd Jsw, dJsw;              //摆动腿的雅可比矩阵
    Eigen::Matrix3d F_fe_rot_sw_W, R_fe_rot_sw_W;            //前后摆动腿在世界坐标系下的旋转矩阵
    Eigen::Vector3d F_fe_pos_sw_W, R_fe_pos_sw_W;            //前后摆动腿在世界坐标系下的位置
    Eigen::Vector3d fe_fl_pos_des_W, fe_fr_pos_des_W, fe_rr_pos_des_W, fe_rl_pos_des_W;         //脚的期望位置
    Eigen::Matrix3d fe_fl_rot_des_W, fe_fr_rot_des_W, fe_rr_rot_des_W, fe_rl_rot_des_W;         //脚的期望旋转矩阵
    Eigen::Vector3d fe_fl_pos_cur_W, fe_fr_pos_cur_W, fe_rr_pos_cur_W, fe_rl_pos_cur_W;         //脚的当前位置
    Eigen::Matrix3d fe_fl_rot_cur_W, fe_fr_rot_cur_W, fe_rr_rot_cur_W, fe_rl_rot_cur_W;         //脚的旋转矩阵
    Eigen::VectorXd q, dq, ddq;         //机器人当前广义坐标位置，速度，加速度
    Eigen::VectorXd Fr_ff;              //MPC计算得出的足端反力
    Eigen::VectorXd delta_ddq;          
    Eigen::VectorXd delta_Fr;
    Eigen::VectorXd eigen_xOpt;         //优化 x=[ddp,fr]
    Eigen::VectorXd eigen_ddq_Opt;      //最终计算得出ddp
    Eigen::VectorXd eigen_fr_Opt, eigen_tau_Opt;    //最终计算得出反力，关节力矩
    Eigen::MatrixXd Q1;
    Eigen::MatrixXd Q2;
    Eigen::VectorXd delta_q_final_kin, dq_final_kin, ddq_final_kin, tauJointRes;    //KinWBC得出的delta_q，dp，ddp，关节力矩
    Eigen::Vector3d pCoMDes, pCoMCur;

    PriorityTasks kin_tasks_walk, kin_tasks_stand;
    // void setQini(const Eigen::VectorXd &qIniDes, const Eigen::VectorXd &qIniCur);
    void setQini(const Eigen::VectorXd &qIni);
    void computeTau();
    void dataBusRead(const DataBus &robotState);
    void dataBusWrite(DataBus &robotState);
    void computeDdq(Pin_KinDyn &pinKinDynIn);
private:
    double timeStep{0.001};
    qpOASES::QProblem QP_prob;
    Eigen::MatrixXd Sf; // floating-base dynamics selection matrix
    Eigen::MatrixXd St_qpV1, St_qpV2; // state selection matrix

    qpOASES::int_t nWSR=100, last_nWSR{0};
    qpOASES::real_t cpu_time=0.1, last_cpu_time{0};
    int qpStatus{0};
    int QP_nv;  //18
    int QP_nc;  //26 20+6
    void copy_Eigen_to_real_t(qpOASES::real_t* target, const Eigen::MatrixXd &source, int nRows, int nCols);
    Eigen::MatrixXd J_base, dJ_base, Jcom;
    Eigen::Vector3d base_pos_des, base_pos, base_rpy_des, base_rpy_cur, hip_link_pos;
    Eigen::Matrix3d hip_link_rot, base_rot;
    Eigen::VectorXd swing_F_fe_pos_des_W, swing_R_fe_pos_des_W;
    Eigen::Vector3d stance_F_fe_pos_cur_W, stance_R_fe_pos_cur_W;
    Eigen::Matrix3d stance_F_fe_rot_cur_W, stance_R_fe_rot_cur_W;
    Eigen::Vector3d stanceDesPos_W;
    Eigen::VectorXd des_ddq, des_dq, des_delta_q, des_q;
    Eigen::VectorXd qIniDes, qIniCur;

    static const int QP_nv_des=18;
    static const int QP_nc_des=26;

    qpOASES::real_t qp_H[QP_nv_des*QP_nv_des];
    qpOASES::real_t qp_A[QP_nc_des*QP_nv_des];
    qpOASES::real_t qp_g[QP_nv_des];
    qpOASES::real_t qp_lbA[QP_nc_des];
    qpOASES::real_t qp_ubA[QP_nc_des];
    qpOASES::real_t xOpt_iniGuess[QP_nv_des];
};


