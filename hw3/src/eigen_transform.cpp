#include "ros/ros.h"
#include <Eigen/Dense>
#include <cmath>
const double PI = 3.14159265359;

Eigen::Vector3d rad2deg(Eigen::Vector3d radians)
{
    // Implement your code here
    radians = radians * 180/PI;
    return radians;
}

Eigen::Vector3d deg2rad(Eigen::Vector3d degrees)
{
    // Implement your code here
    degrees = degrees * PI/180;
    return degrees;
}

Eigen::Quaterniond Euler2Quaternion(Eigen::Vector3d euler)
{
    Eigen::Quaterniond Q;
    Q = Eigen::AngleAxisd(euler.z(), Eigen::Vector3d::UnitZ()) *
        Eigen::AngleAxisd(euler.y(), Eigen::Vector3d::UnitY()) *
        Eigen::AngleAxisd(euler.x(), Eigen::Vector3d::UnitX());
    return Q;
}

Eigen::Vector3d Quaternion2Euler(Eigen::Quaterniond Q)
{
    Eigen::Vector3d Euler(0, 0, 0);
    Euler.x() = atan2(2 * (Q.w() * Q.x() + Q.y() * Q.z()), (1 - 2 * (Q.x() * Q.x() + Q.y() * Q.y())));
    // Implement your code here
    Euler.y() = asin(2 * (Q.w() * Q.y() + Q.z() * Q.x()));
    Euler.z() = atan2(2 * (Q.x() * Q.y() + Q.w() * Q.z()), (1 - 2 * (Q.y() * Q.y() + Q.z() * Q.z())));
    return Euler;
}

Eigen::Quaterniond Hemilton(Eigen::Quaterniond q1, Eigen::Quaterniond q2)
{
    
    double S1 = q1.w();
    double S2 = q2.w();
    Eigen::Vector3d V1 (q1.x(), q1.y(), q1.z());
    Eigen::Vector3d V2 (q2.x(), q2.y(), q2.z());
    Eigen::Vector3d Vpr;
    Vpr = V2*S1 + V1*S2 + V1.cross(V2);
    Eigen::Quaterniond product (S1*S2 - V1.dot(V2), Vpr.x(), Vpr.y(), Vpr.z());
    return product;
}

Eigen::Quaterniond QuadVecRot(Eigen::Quaterniond p, Eigen::Quaterniond q)
{
    //p' = qpq^-1
    Eigen::Quaterniond temp = Hemilton(p,q.conjugate());
    Eigen::Quaterniond Pprime = Hemilton(q,temp);
    return Pprime;
}


int main(int argc, char **argv)
{
    ros::init(argc, argv, "eigen_transform");
    ros::NodeHandle nh;
    ros::Rate loop_rate(1);

    Eigen::Vector3d world_point(1, 0, 0);
    Eigen::Vector3d tmp_point = world_point;
    Eigen::Vector3d euler_angle_deg(0, 0, 90);
    // porpotion of euler_angle
    double s = 0.66667;
    int count = 0;
    std::cout << "Current point position " << world_point.transpose() << std::endl << std::endl;

    while (ros::ok())
    {
        Eigen::Vector3d euler_angle_rad = deg2rad(euler_angle_deg);
        std::cout << "euler_rad: " << euler_angle_rad.transpose() << std::endl;

        Eigen::Quaterniond Q_Total = Euler2Quaternion(euler_angle_rad);
        std::cout << "Quaternion: " <<Q_Total.coeffs().transpose() << std::endl;

        Eigen::Quaterniond Q_tmp = Eigen::Quaterniond::Identity().slerp(s, Q_Total);
        std::cout << "Q_temp: " << Q_tmp.coeffs().transpose() << std::endl;

        Eigen::Vector3d euler_tmp_rad = Quaternion2Euler(Q_tmp);
        Eigen::Vector3d euler_tmp_deg = rad2deg(euler_tmp_rad);
        std::cout << "Apply rotation roll(X): " << euler_tmp_deg.x()
                  << ", pitch(Y): " << euler_tmp_deg.y()
                  << ", yaw(Z): " << euler_tmp_deg.z() << std::endl;

        // Apply rotation to point with quaternion

        // Implement your code here
        Eigen::Quaterniond Q_point(0, world_point.x(), world_point.y(), world_point.z());
        Eigen::Quaterniond Pprime = QuadVecRot(Q_point, Q_tmp);
        std::cout << "Self-Derived world point: " <<Pprime.coeffs().transpose() <<std::endl;

        world_point = Q_tmp * world_point;
        count ++;
        std::cout << "Current point position " << world_point.transpose() << std::endl << std::endl;

        Eigen::Vector3d point_diff = world_point - tmp_point;
        if (point_diff.norm() < 0.1) {
            std::cout << "Current point rotate " << count << " times to origin position!" << std::endl;
            break;
        }
        loop_rate.sleep();
        ros::spinOnce();
    }
    return 0;
}
