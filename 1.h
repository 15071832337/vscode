/*
 * @Author: ‘15071832337’ ‘1418868984@qq.com’
 * @Date: 2022-07-27 10:01:58
 * @LastEditors: ‘15071832337’ ‘1418868984@qq.com’
 * @LastEditTime: 2022-08-04 09:44:32
 * @FilePath: \code(C++)\vscode\1.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
using namespace std;
int jiecheng_jisuan1();
 void a1();
 void a2();
 void a3();
 void a6();
void swap1(int& i,int& j);
void swapshow1();
void swapshow2();
void swap2(int *p,int*q);
void show2();
void swap3(int x,int y);
void a9();
void a7();
void a8();
int jiecheng_jisuan(int x);
void print(const char *cp);
void printshow();
void CompareIntNumberAndIntPointNumber1(int &x,int *y);
void CompareIntNumberAndIntPointNumber2();
void ToChangePoint(int *x,int *y);
int DeGuiDiaoYong(int val);
void DeGuiDiaoYongShow();
void TheValuEqualXiaBiao();
int TheValuEqualXiaBiao(int );
void BeginAndToShowArr();
void MaxSumOfArr();
void ConstUser(int);
void ConstShow();
void MidShow(int,int,int);
void NiXuOutPutArr(int size);
void MoveArr(int);
void SharingCandy(int);
void LightChange();
void TwoArrShow(int, int);
int F(int,int);
void test01();
void test02();

//定义了一个类
class abc
{
private:
    void LightChange();
protected:
    
public:
    int F(int,int);
    void TwoArrShow(int, int);
    void SharingCandy(int);
    int TheValuEqualXiaBiao(int );
};
//继承上类
class bcd:public abc
{
private:
    
public:
    void MidShow(int,int,int);
   
};

class Person1
{

public:
    Person1(int age){
        this->age=age;//this 指向被调用的成员函数所指向的对象 即对象p1
    }//形参的age和属性的age进行了区分
    int age;
    Person1& PersonAddAge(Person1 &p){
        this->age+=p.age;
        return *this;
        
    }
};

class Person2
{
private:
    /* data */
public:
    int m_A;
    int m_B;
    /*
    Person2 operator+(Person2 &p){
        Person2 temp;
        temp.m_A=this->m_A+p.m_A;
        temp.m_B=this->m_B+p.m_B;
        return temp;
    }
    */
};
class Animal
{
private:
    /* data */
public:
    virtual void speak(){
        cout<<"the animal speak"<<endl;
    }
   
};
class Cat:public Animal
{
private:
    /* data */
public:
void speak(){
    cout<<"cat speak"<<endl;
}
    
};
//地址早绑定，在编译阶段就确定了函数的地址
void doSpeak(Animal &animal){
    animal.speak();
}
//为了处理子类和父类的同名函数谁执行的问题，采用
//virtual地址晚绑定










