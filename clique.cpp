#include <iostream>
#include <sstream>
#include <fstream>
#include <queue>
#include "clique.h"
//#include "mpi.h"
mt19937 rand_num;



Net net;
Net clusterindex;
Net clusternodes;
Edge Edge_LIST;
Edge Edge_LIST1;



deque <int> node_pdf_set;
deque <int> PDF_M;

int N_eff=0;
int NN_cut=0;
double N_eff_T = 0;
double power_N =0;

int cliquenumber;
int M_num;

double mean_r1 = 0.0; // S1最大跳跃的位置 time step
double gap1 = 0.0;    //S1最大跳跃
double mean_m1 = 0.0; // S1最大跳跃的位置 average degree

double sq_mean_r1 = 0.0; // S1最大跳跃的平均位置time step
double sq_mean_m1 = 0.0; // S1最大跳跃的平均位置average degree

double sq_gap1 = 0.0;    //S1最大跳跃
double GAP_1 = 0.0;


double gap_real=0;
double r_real=0;

int round_double(double number)
{
    return (number > 0.0) ? floor(number + 0.5) : ceil(number - 0.5);
}


int findroot(int i)
{


    if(ptr[i]<0) {return i;}

    return ptr[i]= findroot(ptr[i]);


}

string int2str(int i)
{
    string s;
    stringstream ss(s);
    ss << i;

    return ss.str();
}

int kCLIQUEFIND(int v1,int v2)
{
    int newnum=0;
    int k1 = net[v1].size();
    int k2 = net[v2].size();
    if (k1 > k2){int temp = v1;v1 = v2;v2 = temp;}

    if (k1>0 && k2>0)
    {
        for(set<int>::iterator it = net[v1].begin();it != net[v1].end(); it++)
        {
            int v3=*it;
            if(net[v2].find(v3) != net[v2].end()){cliquenodes.push_back(v3);newnum++;}
        }
    }
    return newnum;
}
int Com_Clique(int v1,int v2)
{
    int k1=clusterindex[v1].size();
    int k2=clusterindex[v2].size();
    if(k1 && k2)
    {
        if(k1>k2){int temp=v1;v1=v2;v2=temp;}
        int min1=*(clusterindex[v1].begin());
        int max1=*(clusterindex[v1].rbegin());
        int min2=*(clusterindex[v2].begin());
        int max2=*(clusterindex[v2].rbegin());
        if(min1>max2||min2>max1){return -1;}
        for(set<int>::iterator it = clusterindex[v1].begin();it != clusterindex[v1].end(); it++)
        {
            int v3=*it;
            if(v3<min2)continue;
            else if(v3>max2){return -1;}
            if(clusterindex[v2].find(v3) != clusterindex[v2].end()){return v3;break;}
        }
    }
    return -1;
}
void ClusterCombine(int r1,int r2){

    clusternodes[r1].insert(clusternodes[r2].begin(),clusternodes[r2].end());

    clusternodes[r2].clear();
}
int CLUSTER(int newnum,int v1,int v2)
{
    int r1,r2,v3;
    //creat a new 3-clique community
    ptr.push_back(-1);
    r1=cliquenumber;
    //collect the nodes also belonging to other communities
    set<int> tempset;
    for(int i=0;i<newnum+2;i++)
    {
        v3=cliquenodes[i];
        tempset.insert(v3);
    }
    clusternodes.push_back(tempset);
    tempset.clear();
    //if(v1==2 && v2 ==1) cout<<v3<<endl;
    for(int i=1;i<newnum;i++)
    {
        ptr.push_back(r1);
        clusternodes.push_back(tempset);
    }
    //to refresh the community information
    for(int i=0;i<newnum;i++)
    {
        v3=cliquenodes[i];
        //if(v1==2 && v2 ==1) cout<<v3<<endl;
        int comclique=Com_Clique(v1,v3);
        if(comclique!=-1)
        {
            r2=findroot(comclique);
            if(r1!=r2)
            {
                //int ptrr1=ptr[r1], ptrr2=ptr[r2];
                int k1=clusternodes[r1].size();
                int k2=clusternodes[r2].size();
                if(k1>k2)
                {
                    ClusterCombine(r1,r2);
                    ptr[r2]=r1;
                }
                else
                {
                    ClusterCombine(r2,r1);
                    ptr[r1]=r2;
                    r1=r2;
                }

            }
        }
        comclique=Com_Clique(v2,v3);
        if(comclique!=-1)
        {
            r2=findroot(comclique);
            if(r1!=r2)
            {
                int k1=clusternodes[r1].size();
                int k2=clusternodes[r2].size();
                if(k1>k2)
                {
                    ClusterCombine(r1,r2);
                    ptr[r2]=r1;
                }
                else
                {
                    ClusterCombine(r2,r1);
                    ptr[r1]=r2;
                    r1=r2;
                }


            }
        }
        clusterindex[v1].insert(cliquenumber+i);clusterindex[v2].insert(cliquenumber+i);clusterindex[v3].insert(cliquenumber+i);

        if(i>0)
        {
            int number = cliquenumber+i;

            ptr[number] = r1;

        }
    }
    return clusternodes[r1].size();



}
void percolation_32(vector <int> labels_1)
{

    for(int i=0;i<labels_1.size();i++)
    {
        int v1 = labels_1.at(i);
        for(int j=i+1;j<labels_1.size();j++)
        {
            int v2 = labels_1.at(j);
            if(net[v1].find(v2) ==net[v1].end())
            {
                int newnum = kCLIQUEFIND(v1,v2);
                if(newnum>0)
                {
                    cliquenodes.push_back(v1);cliquenodes.push_back(v2);
                    CLUSTER(newnum,v1,v2);
                    cliquenodes.clear();
                }
                cliquenumber += newnum;
                net[v1].insert(v2);
                net[v2].insert(v1);
            }

        }
    }
}




void read()
{

    string finalname;
    ifstream fintemp;
    Edge_LIST1.resize(NN);


    finalname = "./Authors_networks_clique.dat";
    fintemp.open(finalname.c_str());
    for(int i =0;i<NN;i++)
    {
        int x;
        fintemp >> x;
        for(int j=0;j<x;j++)
        {
            int y;
            fintemp >> y;
            Edge_LIST1[i].push_back(y);

        }
    }

}


void shuffle()
{
    Edge_LIST.resize(NC);
    int startpoint = (NN-NC)*rand_num.genrand_real2();
        set <int> mysetnode;
        for(int j=0;j<NC;j++)
    {
         Edge_LIST[j] =  Edge_LIST1[j+startpoint];

        for(int i=0;i<Edge_LIST1[j+startpoint].size();i++)
        {mysetnode.insert(Edge_LIST1[j+startpoint].at(i));
        }

    }
        N_eff = mysetnode.size();
        N_eff_T += N_eff;

}




int percolate_edge_author(FILE *fp1)
{


    double temp_1 = 0.0;
    double gap_1  = 0.0;
    int r1;
    int num_M = 0;

    net.resize(N);
    clusterindex.resize(N);


    int LargestComponentSize=0 ; //the largest (3,2) cluster

    int temps1 = 0;
    cliquenumber=0;





    shuffle();

    for(int i=0;i<NC;++i)
    {

        if(LargestComponentSize <N_eff )
        {

            int label = i;

            int xxx = Edge_LIST[label].size();
            for(int j=0;j<xxx;j++)
            {
                for(int l=j+1;l<xxx;l++)
                {
                    int v1 = Edge_LIST[label].at(j);
                    int v2 = Edge_LIST[label].at(l);

                    if( net[v1].find(v2) ==net[v1].end() && v1!=v2)
                    {

                        int newnum = kCLIQUEFIND(v1,v2);
                        if(newnum>0)
                        {
                            cliquenodes.push_back(v1);cliquenodes.push_back(v2);
                            int LargestComponentSize_temp=CLUSTER(newnum,v1,v2);
                            if(LargestComponentSize_temp >LargestComponentSize)
                            {
                                LargestComponentSize=LargestComponentSize_temp;
                            }
                            cliquenodes.clear();

                        }
                        cliquenumber += newnum;
                        net[v1].insert(v2);
                        net[v2].insert(v1);
                        num_M++;


                    }



                }

            }

        }
        else
        {
            LargestComponentSize=N_eff;

        }


        S1[i] = double (LargestComponentSize);
        M[i] = double (2*num_M)/N_eff;

        gap_1 = double ( LargestComponentSize -temps1);
        temps1 = LargestComponentSize;
        if(gap_1 >temp_1)  {temp_1 = gap_1; r1 = 2*num_M; }
        //if(temp_1+LargestComponentSize >=N_eff || LargestComponentSize >=0.8*N_eff)
        if(temp_1+LargestComponentSize >=N_eff )
        {

            break;
        }

    }
    net.clear();
    clusterindex.clear();
    clusternodes.clear();
    ptr.clear();
    Edge_LIST.clear();


    double qq1 = double (r1) /N_eff;
    setbuf(fp1,NULL);
    fprintf(fp1,"%d  %d  %.1f %.6f\n", N_eff, LargestComponentSize, temp_1, qq1);
    mean_r1 += qq1;

    GAP_1 += temp_1/LargestComponentSize;
    return LargestComponentSize;


}

int main()
{
    rand_num.init_genrand((unsigned)time( NULL ));

    read();
    string outfile = "./result/";
    string str1=int2str((int)NC);

    string str3="distribution_add_clique_";
    string str5=".out";
    string str7="s_edge_";
    string str4=".dat";
    string filename1=outfile+str3+str1+"_"+str5;
    FILE *fp1;
    fp1=fopen(filename1.c_str(),"ab");


    int LargestComponentSize1;
    double ave_LargestComponentSize = 0;

    for(int step=0;step<AVE;step++)
    {

        LargestComponentSize1 = percolate_edge_author(fp1);
        ave_LargestComponentSize += LargestComponentSize1;


        if(step <10)
        {
            string tempstr=int2str((int)step);
            string filename2=outfile+str7+str1+"_"+tempstr+str4;
            FILE *fp2;
            fp2=fopen(filename2.c_str(),"wb");

            for(int i = 0;i<NC;++i)
            {
                if(S1[i]!=0)
                {
                    fprintf(fp2,"%.6f %.6f\n", M[i],S1[i]/LargestComponentSize1); //S1,Gm最大团，spanning cluster
                }
            }
            fclose(fp2);
        }



    }


    FILE *fp3_core;
    fp3_core=fopen("RCL_clique.dat","ab");
    fprintf(fp3_core,"%.6f %.6f %.6f %.6f\n", N_eff_T/AVE,ave_LargestComponentSize/AVE,mean_r1/AVE, GAP_1/AVE); //S1,Gm最大团，spanning cluster
    fclose(fp3_core);

    fclose(fp1);


    return 0;

}

