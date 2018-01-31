#include <iostream>
#include <vector>
#include <numeric>
#include <sys/time.h>
#include "icp.h"
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;



float my_random(void);
Eigen::Matrix3d rotation_matrix(Eigen::Vector3d axis, float theta);
void test_best_fit(void);
void test_icp(void);
void my_random_shuffle(Eigen::MatrixXd &matrix);


unsigned GetTickCount()
{
        struct timeval tv;
        if(gettimeofday(&tv, NULL) != 0)
                return 0;

        return (tv.tv_sec * 1000) + (tv.tv_usec / 1000);
}


int main(int argc, char *argv[]){

    // srand((int)time(0));

    test_best_fit();
    test_icp();


    return 0;
}



///////////////////////////
//  help function

// 0-1 float variables
float my_random(void){
    float tmp = rand()%100;
    return tmp/1000;
}

void my_random_shuffle(Eigen::MatrixXd &matrix){
    int row = matrix.rows();
    vector<Eigen::Vector3d> temp;
    for(int jj=0; jj < row; jj++){
        temp.push_back(matrix.block<1,3>(jj,0));
    }
    random_shuffle(temp.begin(),temp.end());
    for(int jj=0; jj < row; jj++){
        matrix.block<1,3>(jj,0) = temp[jj].transpose();
        // cout << temp[jj].transpose() << endl;
        // cout << "row  " << row << endl;
    }
}


Eigen::Matrix3d rotation_matrix(Eigen::Vector3d axis, float theta){
    axis = axis / sqrt(axis.transpose()*axis);
    float a = cos(theta/2);
    Eigen::Vector3d temp = -axis*sin(theta/2);
    float b,c,d;
    b = temp(0);
    c = temp(1);
    d = temp(2);
    Eigen::Matrix3d R;
    R << a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c),
        2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b),
        2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c;

    return R;
}


void test_best_fit(void){
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(N_pt,3);
    Eigen::MatrixXd B;
    Eigen::MatrixXd C;
    Eigen::Vector3d t;
    Eigen::Matrix3d R;
    Eigen::Matrix4d T;
    Eigen::Vector3d t1;
    Eigen::Matrix3d R1;

    float total_time = 0;
    unsigned start, end;
    float interval;

    for (int i=0; i < N_tests; i++){
        B = A;
        t = Eigen::Vector3d::Random()*translation;

        for( int jj =0; jj< N_pt; jj++){
            B.block<1,3>(jj,0) = B.block<1,3>(jj,0) + t.transpose();
        }

        R = rotation_matrix(Eigen::Vector3d::Random() ,my_random()*rotation);
        B = (R * B.transpose()).transpose();

        B += Eigen::MatrixXd::Random(N_pt,3) * noise_sigma;

        start = GetTickCount();
        T = best_fit_transform(B,A);
        end = GetTickCount();
        interval = float((end - start))/1000;
        total_time += interval;

        C = Eigen::MatrixXd::Ones(N_pt,4);
        C.block<N_pt,3>(0,0) = B;

        C = (T * C.transpose()).transpose();
        t1 = T.block<3,1>(0,3);
        R1 = T.block<3,3>(0,0);

        if(i == 3){
            cout << "position error" << endl << C.block<N_pt,3>(0,0) - A << endl << endl;
            cout << "trans error" << endl << -t1 - t << endl << endl;
            cout << "R error" << endl << R1.inverse() - R << endl << endl;
        }
    }
    cout << "best fit time: " << total_time/N_tests << endl;
}


void test_icp(void){
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(N_pt,3);
    Eigen::MatrixXd B;
    Eigen::MatrixXd C;
    Eigen::Vector3d t;
    Eigen::Matrix3d R;
    Eigen::Matrix4d T;
    Eigen::Vector3d t1;
    Eigen::Matrix3d R1;
    ICP_OUT icp_result;
    std::vector<float> dist;
    int iter;
    float mean;

    float total_time = 0;
    unsigned start, end;
    float interval;

    for (int i=0; i < N_tests; i++){
        B = A;
        t = Eigen::Vector3d::Random()*translation;

        for( int jj =0; jj< N_pt; jj++){
            B.block<1,3>(jj,0) = B.block<1,3>(jj,0) + t.transpose();
        }

        R = rotation_matrix(Eigen::Vector3d::Random() ,my_random()*rotation);
        B = (R * B.transpose()).transpose();

        B += Eigen::MatrixXd::Random(N_pt,3) * noise_sigma;

        // shuffle
        my_random_shuffle(B);

        start = GetTickCount();
        icp_result = icp(B, A, 20,  0.000001);
        end = GetTickCount();
        interval = float((end - start))/1000;
        // cout << "interval" << interval << endl;
        total_time += interval;

        T = icp_result.trans;
        dist = icp_result.distances;
        iter = icp_result.iter;
        mean = std::accumulate(dist.begin(),dist.end(),0.0)/dist.size();

        C = Eigen::MatrixXd::Ones(N_pt,4);
        C.block<N_pt,3>(0,0) = B;

        C = (T * C.transpose()).transpose();
        t1 = T.block<3,1>(0,3);
        R1 = T.block<3,3>(0,0);


        if(i == 3){
            cout << "mean error is " << mean - 6*noise_sigma << endl << endl;
            cout << "icp trans error" << endl << -t1 - t << endl << endl;
            cout << "icp R error " << endl << R1.inverse() - R << endl << endl;
        }

    }
    cout << "icp time: " << total_time/N_tests << endl;

}



