#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<time.h>
using namespace std;

//定义圆周率
#define PI 3.141592657

//定义常数项
#define eplison 1e-3

//定义域
#define X1 -1.000000
#define X2 2.000000

//定义初始种群中有830个个体
#define num 830
//精度需要小数点后六位，编码长度至少需要22位
#define lenOfDNA 22

//繁殖率
#define reproduce_rate 0.7

//变异率
#define mutation_rate 0.25

//定义函数
double* F(double* X){
    int m = _msize(X)/sizeof(double);           //计算坐标个数
    double* Y = new double[m];
    for(int i = 0;i < m;i++){
        Y[i] = X[i]*sin(10*PI*X[i]) +1.0;
    }
    return Y;
}

//随机初始化种群并对其编码
int** initPopulation(){
    int** pop = new int*[num];
    for(int i=0;i < num;i++){
        pop[i] = new int[lenOfDNA];
        for(int j=0;j < lenOfDNA;j++){
            pop[i][j] = rand()%2;              //初始化DNA，注意，每次初始化都是相同的结果，因为用的rand()
        }
    }
    return pop;
}

//解码操作，将DNA转为坐标x
double* decode(int** pop){
    int m = _msize(pop)/sizeof(*pop);                 //计算种群的个体数量，注意_msize是widows独有
    double* X = new double[m];                        //申请数组，代表每个个体的坐标
    for(int i=0;i < m;i++){
        int dna = 0;
        for(int j=0;j < lenOfDNA;j++){
            dna += pop[i][j]*pow(2,j)*1.0;           //对二进制编码加权求和转为十进制
        }
        double x = X1 + dna*(X2-X1)/(pow(2,22)-1.);  //转为对应区间的实数坐标
        X[i] = x;
    }
    return X;
}

//找出最小值,双指针法
double min(double* Y){
    int m = _msize(Y)/sizeof(double);
    double min = Y[0];
    for(int i = 1;i < m;i++){
        min = min<Y[i]?min:Y[i];
    }
    return min;
}

//定义适应度函数,群体中的最大值减去最小值加上一个非常小的数来评估
double* fitness(int** pop){
    double* X = decode(pop);                       //先解码获得坐标X的矩阵
    double* Y = F(X);                               //解坐标对应的Y值
    int m = _msize(Y)/sizeof(double);
    double mini = min(Y);
    double* fit = new double[m];
    for(int i = 0;i < m;i++){
        fit[i] = Y[i] - mini + eplison;
    }
    delete [] X;            //释放空间
    X = NULL;
    return fit;
}

//计算累计概率
double* prop(double* fit){
    int m = _msize(fit)/sizeof(double);
    double sum = 0;
    for(int i = 0;i < m;i++){
        sum += fit[i];
    }
    double* fitProb = new double[m];
    fitProb[0] = fit[0]/sum;
    for(int i = 1;i < m;i++){
        fitProb[i] = fitProb[i-1] + fit[i]/sum;       //概率累加
    }
    delete [] fit;                  //释放空间
    fit = NULL;
    return fitProb;
}

//产生0~1之间的随机数
double getRand(){
    return rand()/(double(RAND_MAX));
}

//二分法查找
int binarySearch(double sr,double* fitProb){
    if(sr<=fitProb[0]){
        return 0;
    }
    int m = _msize(fitProb)/sizeof(double);
    int left = 0;
    int right = m-1;
    int mid = (left+right)/2;
    while(left+1<right){
        if(sr>fitProb[mid]){
            left = mid;
        }else if (sr<fitProb[mid]){
            right = mid;
        }else{
            return mid;
        }
        mid = (left+right)/2;
    }
    return right;
}

//对数组排序
void sort(int left,int right,int* index){
    if(left>=right){
        return;
    }
    int i = left;
    int j = right;
    int base = index[left];
    int temp;
    while(i<j){
        while(index[j]>=base&&i<j){
            j--;
        }
        while(index[i]<=base&&i<j){
            i++;
        }
        if(i<j){
            temp = index[i];
            index[i] = index[j];
            index[j] = temp;
        }
    }
    index[left] = index[i];
    index[i] = base;
    sort(left,i-1,index);
    sort(i+1,right,index);
}

//删除重复元素
int* getIndex(int* index){
    int m = _msize(index)/sizeof(int);
    sort(0,m-1,index);
    int i = 1;
    int sum = 0;
    while(i<m){
        if(index[i]!=index[i-1]){
            sum++;
            index[sum] = index[i];
        }
        i++;
    }
    int* newIndex = new int[sum+1];
    for(int i = 0;i<=sum;i++){
        newIndex[i] = index[i];
    }
    delete [] index;            //释放空间
    return newIndex;
}

//选择函数，适应度高的更容易被选择
int** select(int** pop,double* fitProb){
    int m = _msize(pop)/sizeof(*pop);
    int* sIndex = new int[m];
    for(int i = 0;i < m;i++){
        double sr = getRand();              //产生随机因子
        int index = binarySearch(sr,fitProb);   //选择个体
        sIndex[i] = index;
    }
    int* index = getIndex(sIndex);
    int n = _msize(index)/sizeof(int);
    int** newPop = new int*[n];
    for(int i = 0;i < n;i++){
        newPop[i] = new int[lenOfDNA];
        for(int j = 0;j < lenOfDNA;j++){
            newPop[i][j] = pop[index[i]][j];        //深拷贝
        }
    }
     for(int i = 0;i < m;i++){
        delete [] pop[i];
        pop[i] = NULL;
    }
    delete [] pop;              //释放空间
    pop = NULL;
    delete [] fitProb;          //释放空间
    fitProb = NULL;
    return newPop; 
}

//变异函数
void mutation(int* child){
    if(getRand()<mutation_rate){
        int mutation_point = rand()%lenOfDNA;         //变异点
        child[mutation_point] ^= 1;
    }
}

//繁殖，交叉函数
int** reproduce(int** pop){
    int m = _msize(pop)/sizeof(*pop);
    int** nextPop = new int*[m];
    for(int i = 0;i < m;i++){
        nextPop[i] = new int[lenOfDNA];
    }
    int* child = new int[lenOfDNA];
    int* mother = new int[lenOfDNA];
    for(int i = 0;i < m;i++){
        for(int j = 0;j < lenOfDNA;j++){
            child[j] = pop[i][j];
        }
        if(getRand()<reproduce_rate){
            int other = rand()%m;
            for(int j = 0;j < lenOfDNA;j++){
                mother[j] = pop[other][j];
            }
            int cross_point = rand()%lenOfDNA;
            for(int j = cross_point;j < lenOfDNA;j++){
                child[j] = mother[j];
            }
        }
        mutation(child);
        for(int j = 0;j < lenOfDNA;j++){
            nextPop[i][j] = child[j];
        }
    }
    delete [] child;
    child = NULL;
    delete [] mother;
    mother = NULL;
    for(int i = 0;i < m;i++){
        delete [] pop[i];
        pop[i] = NULL;
    }
    delete [] pop;
    pop = NULL;
    return nextPop;
}

//迭代N次
void iteration(int N){
    int** pop = initPopulation();           //初始化种群
    for(int i = 0;i < N;i++){
        double* fit = fitness(pop);         //适应度评价
        double* fitProb = prop(fit);        //获得累积概率
        int** newpop = select(pop,fitProb); //选择
        pop = reproduce(newpop);  //交叉变异
    }
    cout<<"last iteration:"<<endl;
    double* X = decode(pop);
    double* Y = F(X);
    int m = _msize(X)/sizeof(double);
    for(int i = 0;i < m;i++){
        cout<<"x:"<<X[i]<<endl;
    }
    int n = _msize(Y)/sizeof(double);
    for(int i = 0;i < m;i++){
        cout<<"y:"<<Y[i]<<endl;
    }
}

int main(){
    iteration(50);
    getchar();
    return 0;
}