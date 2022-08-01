/*
 * @Author: ‘15071832337’ ‘1418868984@qq.com’
 * @Date: 2022-07-27 10:01:58
 * @LastEditors: ‘15071832337’ ‘1418868984@qq.com’
 * @LastEditTime: 2022-08-01 11:07:35
 * @FilePath: \code(C++)\vscode\1.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
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




