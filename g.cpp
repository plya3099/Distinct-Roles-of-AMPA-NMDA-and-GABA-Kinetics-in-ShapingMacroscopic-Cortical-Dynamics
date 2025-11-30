//
// Created by apjas on 2024/4/12.
//
#include <iostream>
#include <fstream>
#include "math.h"
//#include "iostream.h"
#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "time.h"
#include <string.h>
//#include<ctime>
#include <random>
using namespace  std;

#define pi 3.14159265358979303846

int n = 12500;
int ne = 10000;
int ni = 2500;
double epsilon = 1e-6;


namespace Random {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<double> dis(0.0, 1.0);
}



struct Node{
    int label;  //神经元类型，1为兴奋性神经元，0为抑制性神经元
    int* same_list; // 与自身相同的神经元类型相连接节点的列表
    int* diff_list;//  与自身类型不相同的神经元类型相连接的列表
    int  k;//度
    int k2;
    //double *g1,*g2; //突触强度
    double byq; //进入不应期时间
    double spike_times;
    double v;//电压
    bool firing ;
    int isi_count;    // ISI数量
    double isi_mean;  // ISI均值
    double isi_m2;    // ISI方差中间量（Welford算法）

    //bool block;
}neuron;

void init_stats(struct Node *stats) {
    for (int i = 0; i < n; i++) {
        stats[i].byq = -1.0; // 表示尚未发放
        stats[i].isi_count = 0;
        stats[i].isi_mean = 0.0;
        stats[i].isi_m2 = 0.0;
    }
}




double random_float() {
    return Random::dis(Random::gen);
}

double getGaussianNoise(double mean, double stddev) {
    std::normal_distribution<double> dist(mean, stddev);
    return dist(Random::gen);
}


void update_stats(struct Node  *stats, int neuron_id, double current_time) {
    if (stats[neuron_id].byq < 0.0||current_time<100.0) {
        // 第一次发放，只记录时间
        stats[neuron_id].byq = current_time;
        return;
    }

    // 计算新的ISI
    double isi = current_time - stats[neuron_id].byq;
    stats[neuron_id].byq = current_time;

    // 更新ISI统计量（Welford算法）
    stats[neuron_id].isi_count++;
    double delta = isi - stats[neuron_id].isi_mean;
    stats[neuron_id].isi_mean += delta / stats[neuron_id].isi_count;
    double delta2 = isi - stats[neuron_id].isi_mean;
    stats[neuron_id].isi_m2 += delta * delta2;
}


// 计算网络平均CV
double calculate_network_cv(struct Node *stats) {
    double total_cv = 0.0;
    int valid_neurons = 0;

    for (int i = 0; i < n; i++) {
        if (stats[i].isi_count < 1) continue; // 至少需要2次脉冲（1个ISI）

        // 计算标准差
        double variance = stats[i].isi_m2 / stats[i].isi_count;
        double std_dev = sqrt(variance);

        // 计算CV
        double cv = stats[i].isi_mean > 0.0 ? std_dev / stats[i].isi_mean : 0.0;
        if (cv > 0.0) {
            total_cv += cv;
            valid_neurons++;
        }
    }

    return valid_neurons > 0 ? total_cv / valid_neurons : 0.0;
}




int cauchy_random(int x0, double gamma) {
    
    double u = Random::dis(Random::gen);

    double result = x0 + gamma * std::tan(pi * (u - 0.5));
    if (result < 0.0) result = -result;
    return static_cast<int>(std::round(result));
}



int random_num(int min,int max) {
    if (min < 0) {
        return ceil(random_float() * max);
    } else {

        return min+floor(random_float() * (max - min));
    }

}

inline double B(double v) {

    return 1/(1+0.280112*exp(-0.062*v));

}


int main(int argc,char **argv)
{
  long kin_sum=0;

    char *pt;
    ///double control_tau = 0.0;
    double tau_ampa = 0.2;
    double tau_nmda =8;
    double tau_gaba = 0.25;
    double Jset=1;
    double Sapma=0;
    double Snmda=0;
    double Sgaba=0;
    double i00=1;
    double g00=1;

    //double i_ext[10000];
/*

    for(int i=0;i<argc;i++){
        if(strstr(argv[i],"-tau_ampa=")==argv[i]){
            pt = strstr(argv[i],"=")+1;
            tau_ampa = atof(pt);
            tau_ampa = tau_ampa/20.0;
        }
        if(strstr(argv[i],"-tau_nmda=")==argv[i]){
            pt = strstr(argv[i],"=")+1;
             tau_nmda  = atof(pt);
            tau_nmda = tau_nmda/20.0;   
        }
	 if(strstr(argv[i],"-tau_gaba=")==argv[i]){
            pt = strstr(argv[i],"=")+1;
            tau_gaba  = atof(pt);
            tau_gaba = tau_gaba/20.0;
        }



	  if(strstr(argv[i],"-i=")==argv[i]){
            pt = strstr(argv[i],"=")+1;
            i00= atof(pt);
            
        }
	if(strstr(argv[i],"-j=")==argv[i]){
		pt = strstr(argv[i],"=")+1;
		Jset=atof(pt);	
	}

        if(strstr(argv[i],"-g=")==argv[i]){
                pt = strstr(argv[i],"=")+1;
                g00=atof(pt);
        }



	        if(strstr(argv[i],"-sapma=")==argv[i]){
                pt = strstr(argv[i],"=")+1;
                Sapma=atof(pt);
                }
        if(strstr(argv[i],"-snmda=")==argv[i]){
                pt = strstr(argv[i],"=")+1;
                Snmda=atof(pt);
                }

        if(strstr(argv[i],"-sgaba=")==argv[i]){
                pt = strstr(argv[i],"=")+1;
                Sgaba=atof(pt);
                }
    }

*/
cout<<tau_ampa<<" , "<<tau_nmda<<" , "<<i00<<endl;
    int k0=1000;//平均度
    double v_p=70,v_r=-70;
// The parameter about the input current:
    double f0=1.0,big_jee=0.275*f0,big_jii=0.953939*f0,big_jie=0.3*f0,big_jei=0.96286*f0;
    double i_0e=i00*g00,i_0i=i00;//i_0i=i_0e/1.02;
    double tau=0.01,delta=0.0;
    double t_rv=2.0/v_p;
// a sliding time window of size
    double delta_0ii,delta_0ee;
    double i0e,i0i,delta_kii,delta_kee,j0ii,j0ee,j0ie,j0ei,big_jie0,big_jei0;
    double Time =600.0;
    double twin = 0.0;
    double dt=0.0001;
    double delta_t = 0.005;

    double v_all= 0.0;
    double v_e=0.0,v_i=0.0;
    double v_aver =0.0;
    double r_all=0.0;
    double r_e=0.0;
    double r_i=0.0;
    double r_aver=0.0;



    double compute_times=0.0;
    double normal_num = 0.0;
    double normal_e = 0.0;
    double normal_i = 0.0;


    struct Node neurons[n];

    delta_0ee=0.1/f0;
    delta_0ii=0.3/f0;
    i0e=0.1*sqrt(1000);//i_0e*sqrt(k0);
    i0i=i0e/1.02;//1.02;//i_0i*sqrt(k0);
    delta_kee=delta_0ee*sqrt(k0);
    delta_kii=delta_0ii*sqrt(k0);
    big_jie0=big_jie*k0/k0;
    big_jei0=big_jei*k0/k0;
    j0ee=big_jee*sqrt(k0);
    j0ii=-big_jii*sqrt(k0);
    j0ie=big_jie0*sqrt(k0);
    j0ei=-big_jei0*sqrt(k0);



    // bool adjacency_matrix[n][n]={false}; //连接矩阵
    bool nodes_list[n];
    for(int j=0;j<n;j++){nodes_list[j]=false;}

    char file_name[60];

   // FILE *fp_e;
    //FILE *fp_i;

    sprintf(file_name, "data.csv");
    ofstream  fp_data(file_name,ios::out);

  // sprintf(file_name, "single.csv");
   // ofstream  fp_single(file_name,ios::out);

   // sprintf(file_name, "single2.csv");
  //  ofstream  fp_single2(file_name,ios::out);

  //  fp_single<<"v,i_ampa,i_nmda,i_gaba,s_ampa,s_nmda,s_gaba"<<endl;
   // fp_single2<<"v,i_ampa,i_nmda,i_gaba,s_ampa,s_nmda,s_gaba"<<endl;
    fp_data<<"t,rall,vall,re,ri,ve,vi,saee,snee,sgei,saie,snie,sgii"<<endl;
    if(!fp_data.is_open()){
        cout<<"文件打开失败"<<endl;
        return  0;
    }
    char ch[50];
   sprintf(ch, "data_e.csv");

    ofstream fp_e(ch,ios::out);
  sprintf(ch, "data_i.csv");

    ofstream fp_i(ch,ios::out);



    /*

    sprintf(file_name, "data_e：ampa:%.4f,nmda:%.4f.txt",tau_ampa,tau_nmda);
    fp_e = fopen(file_name, "w");


    sprintf(file_name, "data_i.txt");
    fp_i = fopen(file_name, "w");

  */
   // sprintf(file_name, "a-%ld-data%.6f.txt",time(0),control_tau);
 //   vr=fopen(file_name, "w");
    // int count = 0;


    int kee=0,kii=0;


    /********初始化神经网络******/
    // int num=0;
    int node;
    srand(time(NULL));//重置随机数种子
double rnum;


    long allK=0;

    for(int i=0;i<n;i++) {

      

        //
        /*******初始化连接******/

        /**初始化相同类型的神经元之间的连接**/
        if (i < ne) {  //初始化神经元的度，兴奋性神经元和抑制性神经元的度分布的尺度参 数不一样
            neurons[i].label = 1;

            neurons[i].k =  cauchy_random(k0, 10);
            while (neurons[i].k >= 10000){neurons[i].k =  cauchy_random(k0, 10);}
           // cout<<neurons[i].k <<endl;
          neurons[i].k2 =250;

           // cout<<i<<":"<<neurons[i].k2<<endl;
           // i_ext[i] = cauchy_random2(i0e,0.0003);
            //要确保生成的度在设置的数目之间
           // while (neurons[i].k < 1 || neurons[i].k > ne) { neurons[i].k = cauchy_random(k0, delta_kee); } //
            //记录该节点与相同类型的神经元相连的度
        } else {
            //i_ext[i] = cauchy_random2(i0i,0.0003);
          neurons[i].k =  cauchy_random(k0, 10);
          neurons[i].k2 =  4000;
            if (neurons[i].k >= 2500){neurons[i].k =  cauchy_random(k0, 10);}
        }


        if (i < ne) {
            //e to e
            neurons[i].same_list = new int[neurons[i].k+10];
            int count=0;
            while (count<neurons[i].k) {
                node = random_num(-1,ne);
                while (nodes_list[node]!=false||node==i) {
                    node = random_num(-1,ne);
                }
                    neurons[i].same_list[count]=node;
                nodes_list[node] = true;
                if (node==11000){kee++;}
                allK +=1;
                    count++;

            }
            neurons[i].k = count;



        } else {
            // i to i
            neurons[i].same_list = new int[neurons[i].k+100];
            int count=0;
            while (count<neurons[i].k) {
                node = random_num(ne,n);
                while (nodes_list[node]||node==i) {
                    node = random_num(ne,n);

                }

                nodes_list[node] = true;
                neurons[i].same_list[count]=node;
                if (node==12000){kii++;}
                allK +=1;
                count++;

            }

            neurons[i].k = count;
            kin_sum += neurons[i].k;

        }




        if (  i < ne) {
            //e to i

            neurons[i].diff_list = new int[neurons[i].k2+100];
            int count=0;
            while (count<neurons[i].k2) {

                node = random_num(ne,n);
                while (nodes_list[node]!=false||node==i) {
                    node = random_num(ne,n);
                }
                nodes_list[node] = true;
                neurons[i].diff_list[count]=node;
                allK +=1;

                if (node==12000){kee++;}
                count++;


            }
            neurons[i].k2 = count;
            //cout<<neurons[i].k2<<endl;

        } else {
            //i to e
            neurons[i].diff_list = new int[neurons[i].k2+100];
            int count=0;
            while (count<neurons[i].k2) {
                node = random_num(-1,ne);
                while (nodes_list[node]!=false||node==i) {
                    node = random_num(-1,ne);
                }
                nodes_list[node] = true;
                neurons[i].diff_list[count]=node;
                allK +=1;
                count++;

            }
            neurons[i].k2 = count;

        }

//	printf("passx5\n");
        /**初始化网络的初始电压**/
        neurons[i].v = 10+20*( (double) ((double) rand() / (double) RAND_MAX));

        neurons[i].byq = 0.0;
        neurons[i].spike_times = 0.0;


        //printf("g init finished\n");

        neurons[i].firing = false;
        // neurons[i].block = false;

        for (int j = 0; j < n; j++) { nodes_list[j] = false; }
	//cout<<i<<" : "<<neurons[i].k<<" , "<<neurons[i].k2<<endl;

    }
//    printf("%d\n",
// cout<<allK/25000<<endl;
    printf("time evolution2025\n");
  cout<< kin_sum<<endl;
  //  cout<<"no.5000 Kin from e:"<<kee<<" , Kin from i : "<<kii<<endl;

cout<<i_0e<<","<<i_0i<<endl;
    //bool bock = false;
    //int point = 0;

    int i;
    double I_ampa[14000]={0.0};
    double I_nmda[14000]={0.0};
    double I_gaba[14000]={0.0};
   // double I_in[5000]={0.0};
    double Srec_ampa[14000]={Sapma};
    double Srec_nmda[14000]={Snmda};
    double Srec_gaba[14000]={Sgaba};
    //double S_ampa[5000]={0.0};
   // double S_nmda[5000]={0.0};
    //double S_gaba[5000]={0.0};
    double  jaee=0.6;// 0.6;//0.01897366;
    double  jnee=0.178;//0.178;//0.00562289;//0.005622;
    double  jgei=0.0155;;//0.076;//0.011384196;//0.0004905;
    ;//0.003964034;
    double  jaie=0.597014925373134; //0.597014925373134;//0.01896;//0.005622;
    double  jnie= 0.177114427860697;// 0.177114427860697;
    double  jgii= 0.0153465346534653;//0.0759240759240759;//0.011384;//0.000464034;



    //double g_rec_ampa=tau_ampa,g_nmda=,g_gaba=1.3;

   // double tau_ampa = 0.1;
    //double tau_nmda = 5.0;
    //double tau_gaba = 0.1;

    cout<<"ee"<<j0ii/k0<<endl;
    /*****时间演化******/

    init_stats(neurons);

    bool noise = false;
    for(double t=0;t<Time;t=t+dt){

/*
        if (abs(t - 250) < epsilon) {

           //i00=0.1;
            cout<<"p1111"<<endl;
            //i_0e=i00;
           // i_0i=i00;
           // i0e=i_0e*sqrt(k0);
          //  i0i=i_0i*sqrt(k0);
            noise = true;
           // tau_gaba=0.315;
            //tau_nmda = 10;
        }

        if (abs(t - 550) < epsilon) {
            cout<<"p222"<<endl;
            //i00=0.01;
           // cout<<"p2222"<<endl;
          //  i_0e=i00;
          //  i_0i=i00;
            //i0e=i_0e*sqrt(k0);
           noise = false;
           // tau_gaba=0.325;
          //  i0i=i_0i*sqrt(k0);
            //tau_nmda = 3.5;
        }
*/


        if (t - twin > delta_t) {
            twin = t;

            if (t >0.0) {

                /*开始统计网络信息*/
                /***整体网络信息统计***/
                //计算平均发放率


                //int num = 0;
                v_all = 0;
                v_e=0.0;
                v_i=0.0;
                r_e=0.0;
                r_i=0.0;
                r_all = 0;
                compute_times += 1;
                normal_num = 0;
                normal_e = 0;
                normal_i = 0;
                double s_n_ee,s_a_ee,s_g_ei;
                double s_n_ie,s_a_ie,s_g_ii;
                //  compute_times +=1;
                for (i = 0; i < ne; i++) {
                    /*
                    if(neurons[i].byq>(2.0/v_p)){
                        v_all = v_all +  neurons[i].v;
                        num = num +1;
                    }
                    */
                    r_all = r_all + +neurons[i].spike_times / (n * delta_t);
                    r_e = r_e + neurons[i].spike_times / (ne * delta_t);
                    neurons[i].spike_times = 0;
                    if ((t - neurons[i].byq) > 0.2) {
                        v_all = v_all + neurons[i].v;
                        v_e = v_e + neurons[i].v;
                        normal_num += 1;
                        normal_e += 1;
                    }
                    s_n_ee += I_nmda[i];
                    s_a_ee += I_ampa[i];
                    s_g_ei += I_gaba[i];

                }
                s_n_ee = s_n_ee /5000;
                s_a_ee = s_a_ee/5000;
                s_g_ei = s_g_ei/5000;

                for (i = ne; i < n; i++) {
                    /*
                    if(neurons[i].byq>(2.0/v_p)){
                        v_all = v_all +  neurons[i].v;
                        num = num +1;
                    }
                    */
                    r_all = r_all + +neurons[i].spike_times / (n * delta_t);
                    r_i = r_i + neurons[i].spike_times / (ni * delta_t);
                    neurons[i].spike_times = 0;
                    if ((t - neurons[i].byq) > 0.2) {
                        v_all = v_all + neurons[i].v;
                        v_i = v_i + neurons[i].v;
                        normal_num += 1;
                        normal_i += 1;
                    }
                    s_n_ie += I_ampa[i];
                    s_a_ie += I_nmda[i];
                    s_g_ii += I_gaba[i];


                }
                s_n_ie = s_n_ie /5000;
                s_a_ie = s_a_ie/5000;
                s_g_ii = s_g_ii/5000;



                if (normal_num > 0) { v_all = v_all / normal_num; } else { v_all = 0; }
                if (normal_e > 0) { v_e = v_e / normal_e; } else { v_e = 0; }
                if (normal_i > 0) { v_i = v_i / normal_i; } else { v_i = 0; }
                //v_aver = ((compute_times - 1) * v_aver + v_all) / compute_times;
                //r_aver = ((compute_times - 1) * r_aver + r_all) / compute_times;




                fp_data<<t<<","<<50*r_all<<","<<v_all <<","<<50*r_e<<","<<50*r_i<<","<<v_e<<","<<v_i<<",";
                fp_data<<s_a_ee<<","<<s_n_ee<<","<<s_g_ei<<","<<s_a_ie<<","<<s_n_ie<<","<<s_g_ii<<endl;
                //r_aver = r_aver / compute_times;



            }else{
                for (i = 0; i < n; i++) {
                    neurons[i].spike_times = 0;
                }
            }

        }

        /***遍历所有神经元***/


     //cout<<"p1"<<endl;
        /***遍历所有神经元***/
        for(i=0;i<ne;i++){
            /**外部驱动电流的计算**/


            I_nmda[i] =jnee*(neurons[i].v)*Srec_nmda[i]/sqrt(k0);//0.00562289
            I_ampa[i] =jaee*(neurons[i].v)*Srec_ampa[i]/sqrt(k0);//0.01897366
            I_gaba[i] =jgei*(neurons[i].v+70)*Srec_gaba[i]/sqrt(k0);//0.0004905

            if(t-neurons[i].byq-0.2>=0.0) {

                    neurons[i].v=neurons[i].v+(neurons[i].v*neurons[i].v - I_nmda[i] - I_gaba[i]  - I_ampa[i] + i0e  )*dt;

            }
            Srec_ampa[i] = Srec_ampa[i]  + (-Srec_ampa[i]/tau_ampa )*dt;
            Srec_nmda[i] = Srec_nmda[i]  + (-Srec_nmda[i]/tau_nmda )*dt ;
            Srec_gaba[i] = Srec_gaba[i] + (-Srec_gaba[i]/tau_gaba)*dt ;

        }

	   for(i=ne;i<n;i++){
            /**外部驱动电流的计算**/

           I_nmda[i] = jnie*(neurons[i].v)*Srec_nmda[i]/sqrt(k0);//0.005622
           I_ampa[i] = jaie*(neurons[i].v)*Srec_ampa[i]/sqrt(k0); // 0.01896
           I_gaba[i] =jgii*(neurons[i].v+70)*Srec_gaba[i]/sqrt(k0);//0.000464034

            if(t-neurons[i].byq- 0.2>=0.0)
            {


                    neurons[i].v=neurons[i].v+(neurons[i].v*neurons[i].v - I_nmda[i] - I_gaba[i]  - I_ampa[i] + i0i )*dt;


            }

           Srec_ampa[i] = Srec_ampa[i]  + (-Srec_ampa[i]/tau_ampa )*dt;
           Srec_nmda[i] = Srec_nmda[i]  + (-Srec_nmda[i]/tau_nmda )*dt ;
           Srec_gaba[i] = Srec_gaba[i] + (-Srec_gaba[i]/tau_gaba)*dt ;

        }



        //cout<<"p2"<<endl;
        for(i=0;i<n;i++)
        {
            if(neurons[i].v>=v_p)
            {
                /*该神经元处于激发状态*/
                update_stats(neurons, i, t);
                neurons[i].v = v_r;
                //neurons[i].byq = t;
                neurons[i].firing = true;
                neurons[i].spike_times = neurons[i].spike_times +1;
                if(i<ne){
                   // S_ampa[i] = S_ampa[i]  + (-S_ampa[i]/tau_ampa )*dt + 1;
                    //S_nmda[i] = S_nmda[i]  + (-S_nmda[i]/tau_nmda )*dt + 1;
                    //S_gaba[i] = S_gaba[i] + (-S_gaba[i]/tau_gaba)*dt ;
                   // fprintf(fp_e,"%d ",i);
                  fp_e<<i<<",";
                }else{
                    //S_ampa[i] = S_ampa[i]  + (-S_ampa[i]/tau_ampa )*dt ;
                    //S_nmda[i] = S_nmda[i]  + (-S_nmda[i]/tau_nmda )*dt ;
                   // S_gaba[i] = S_gaba[i] -(S_gaba[i]/tau_gaba)*dt + 1;
                    //fprintf(fp_i,"%d ",i);
                   fp_i<<i<<",";
               }


            }
            //Srec_ampa[i] = 0.0;
            //Srec_nmda[i] = 0.0;
            //Srec_gaba[i] = 0.0;
        }

        //cout<<"p3"<<endl;

        for(i =0;i<n;i++){
            if((t-neurons[i].byq)>( 0.1)&&neurons[i].firing)
            {
                neurons[i].firing = false;

                //neurons[i].byq = 1000000.0;
                if(i<ne){
                    /*不同类型的神经元产生的驱动电流不同*/

                    for(int j=0;j<neurons[i].k;j++){
                        node = neurons[i].same_list[j];
                        Srec_ampa[node] = Srec_ampa[node] + 1; //0.00341526*5;
                        Srec_nmda[node] = Srec_nmda[node] + 1;//0.001124578*5;
                       // neurons[node].v  =neurons[node].v +j0ee/k0;
                    }

                    for (int j = 0; j < neurons[i].k2; j++) {
                        node = neurons[i].diff_list[j];
                        Srec_ampa[node] = Srec_ampa[node] + 1;// 0.003794732*5;
                        Srec_nmda[node] = Srec_nmda[node] + 1;//0.0011384196*5;
                        //neurons[node].v  = neurons[node].v  +  j0ie/k0;

                    }
/*
                    for(int j=0;j<k0;j++){
                        node = neurons[i].diff_list[j];
                        neurons[node].v =neurons[node].v +neurons[i].res*j0ie/k0;

                    }
                    */

                }else{

                    for(int j=0;j<neurons[i].k;j++){
                        node = neurons[i].same_list[j];
                       // neurons[node].v =neurons[node].v  + j0ii/k0;
                       Srec_gaba[node] = Srec_gaba[node] + 1;// 0.0024525*0.2;;
                    }


                    for (int j = 0; j <  neurons[i].k2; j++) {
                        node = neurons[i].diff_list[j];
                        Srec_gaba[node] = Srec_gaba[node] + 1;//0.00232017*0.2;
                       // neurons[node].v  = neurons[node].v + j0ei/k0;
                    }

                   /*
                    for(int j=0;j<k0;j++){
                        node = neurons[i].diff_list[j];
                        neurons[node].v =neurons[node].v +neurons[i].res*j0ei/k0;
                    }
                 */

               }
            }




        }
 //fp_single<<neurons[5000].v<<","<<I_ampa[5000]<<","<<I_nmda[5000]<<","<<I_gaba[5000]<<","<<Srec_ampa[5000]<<","<<Srec_nmda[5000]<<","<<Srec_gaba[5000]<<endl;
 //fp_single2<<neurons[12000].v<<","<<I_ampa[12000]<<","<<I_nmda[12000]<<","<<I_gaba[12000]<<","<<Srec_ampa[12000]<<","<<Srec_nmda[12000]<<","<<Srec_gaba[12000]<<endl;
        fp_e<<endl;
        fp_i<<endl;
       // cout<<"p5"<<endl;
       // fprintf(fp_e, "\n");
        //fprintf(fp_i, "\n");


    }


    sprintf(file_name, "isi.csv");
    ofstream  fp_isi(file_name,ios::out);
    double avg_cv = calculate_network_cv(neurons);
    fp_isi<<avg_cv<<endl;

  fp_e.close();
    fp_i.close();
    fp_data.close();
    //fclose(fp_data);

    return 0;
}
