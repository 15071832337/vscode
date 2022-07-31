/*
 * @Author: ‘15071832337’ ‘1418868984@qq.com’
 * @Date: 2022-07-22 11:39:19
 * @LastEditors: ‘15071832337’ ‘1418868984@qq.com’
 * @LastEditTime: 2022-07-31 16:28:15
 * @FilePath: \vscode\1.cpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include <iostream>
#include "1.h"
#include <cmath>


using namespace std;
const int MAX = 4;

//定义获取数组长度 sizeof会报错

 
int main ()
{

   
   TwoArrShow(3,3);
   return 0;
}

//利用指针偏移进行数组操作
void a1(){
   int  var[MAX] = {10, 100, 200,300};
   int  *ptr;
 
   // 指针中的数组地址
   ptr = var;
   for (int i = 0; i < MAX; i++)
   {
      cout << "var[" << i << "] the address is ";
      cout << ptr << endl;
      cout << "var[" << i << "] the number is ";
      cout << *ptr << endl;
 
      // 移动到下一个位置
      ptr=ptr+1;
   }
}
//利用数组本身进行数组操作

void a2(){
   int var[MAX] = {10,20,30,40};
   for (int i = 0; i < MAX; i++)
   {
      cout << "var[" << i << "] the address is ";
      cout<<&var[i]<<endl;
      cout << "var[" << i << "] the number is ";
      cout<<var[i]<<endl;
   }
   
}
//总的字节数除以一个元素所占的字节数就是数组的长度
void a3(){
   int  var[5] = {10, 100, 200, 300,500};
   cout<<sizeof(var) <<endl;//给出总的字节数
    cout<<sizeof(var[0])<<endl;//给出一个元素的字节数
   cout<<sizeof(var) / sizeof(var[0])<<endl;
}

//字符数组和字符指针的使用，这里利用函数a4(),a5();编写函数，返回字符c在字符串s中的首次出现的位置
/*
*/
void a6(){
   char s[]={'e','r','t','4','\0'};
   cout<<s[0]<<endl;
   //对字符数组s进行值的修改，全部修改成a
   for (int i = 0; i < sizeof(s)/sizeof(s[0])-1; i++)
   {
      s[i]='a';
   }
   cout<<s<<endl;
   //下面利用指针进行上述的修改操作
   //1.we will firestly define a point signed char and  named t 
   char *t;
   //2.we will put the point of s to t and if the *t is '\0',the for end
   for ( t = s; *t!=s[sizeof(s)/sizeof(s[0])-1]; t++)
   {
      *t='b';//3.the char 'b' given *t because the point is moved,so the *t is changed in for
   }
   cout<<s<<endl;
   cout<<sizeof(s)/sizeof(s[0])-1<<endl;
}

//引用变量是一个别名，也就是说，它是某个已存在变量的另一个名字。一旦把引用初始化为某个变量，就可以使用该引用名称或变量名称来指向变量。
//用引用来交换两个数值


void swapshow1(){
   int a = 10;
   int b = 20;
   cout<<"the numbers is no changed:"<<a<<" "<<b<<endl;
   swap(a,b);
   cout<<"the numbers is  changed:"<<a<<" "<<b<<endl;


}
void swap1(int &i,int &j){
   int temp;
   temp=i;
   i=j;
   j=temp;
   cout<<"the numbers is changed :"<<i<<" "<<j<<endl;
}

void swapshow2(){
   int  x=2;
   int y=4;
   swap1(x,y);
}
//利用指针进行两个数的交换操作，其思路是定义两个整型数，定义两个指针，然后进行指针操作
void show2(){
   int a=10;
   int b=20;
   int *p=&a;
   int *q=&b;
   cout<<"the numbers is no changed "<<a<<" "<<b<<endl;
   swap2(p,q);
   cout<<"the numbers  is changed "<<a<<" "<<b<<endl;
   //这里可以如果用引用的话，可以不用定义整型指针 p,q,可以直接把&a,&b传入swap2
}
void swap2(int *p,int *q){
   int temp;
   temp=*p;
   *p=*q;
   *q=temp;
}
//通过使用引用来替代指针，会使 C++ 程序更容易阅读和维护。C++ 函数可以返回一个引用，方式与返回一个指针类似。
//当函数返回一个引用时，则返回一个指向返回值的隐式指针。这样，函数就可以放在赋值语句的左边。例如，请看下面这个简单的程序：

//下面试验一下这种交换数据的操作,该操作并不能交换两个数的数值，它既不是指针进行地址和内容的操作，也不是引用直接进行操作，这种复制情况只能在函数中交换，函数结束时，临时变量消失。

void a9(){
   int x=10;
   int y=20;
   cout<<"the numbers is no changed :"<<x<<" "<<y<<endl;
   swap3(x,y);
   cout<<"the numbers is  changed? :"<<x<<" "<<y<<endl;
}
void swap3(int x,int y){
   int temp;
   temp=x;
   x=y;
   y=temp;
}
//above example,the answer shown is no changed 

//6.4 编写一个与用户交互的函数，要求用户输入一个数字，计算生成该数字的阶乘。
//第一步传入参数为输入参数，定义阶乘函数，返回一个整数结果
int jiecheng_jisuan(int x){
   
   int i;
   int answer=1;
   if (x>=1)
   {
         for ( i = 1; i < x+1; i++)
      {
         answer=i*answer;
      }
      //return answer;
   cout<<"the answer :"<<answer<<endl;
   }
   else{
      cout<<"the input number x is not included the field of x>=1, you should repeatedly input:"<<endl;
   }
   return 0;
}

int jiecheng_jisuan1(){
   cout<<"Please input a int number x(x>=1):"<<endl;
   int x;
   cin>>x;
   int i;
   int answer=1;
   if (x>=1)
   {
         for ( i = 1; i < x+1; i++)
      {
         answer=i*answer;
      }
      //return answer;
   cout<<"the answer :"<<answer<<endl;
   }
   else{
      cout<<"the input number x is not included the field of x>=1, you should repeatedly input:"<<endl;
   }
   return 0;
}

//字符的输出函数
void print(const char *cp){
   if (cp)
   {
      while (*cp)
      {
         cout<<*(cp++);
      }
      
   }
   
}
void printshow(){
   char s[]={'a','r','d','\0'};
   print(s);
}

//编写一个练习，接收两个参数，一个是int的数，一个是int型指针，函数比较int的值和指针所指的值，返回较大的那个。
void CompareIntNumberAndIntPointNumber1(int &x,int *y){
   if (x>*y)
   {
      cout<<"the bigger number is:"<<x<<endl;
   }
   else if (x<*y)
   {
      cout<<"the bigger number is:"<<*y<<endl;
   }
   else{
      cout<<"the two numbers is same"<<endl;
   }
}
//这里为了方便，我将输入的两个整数，一个直接利用引用，另一个把输入内容的地址交给指针
void CompareIntNumberAndIntPointNumber2(){
   int x;
   int y;
   int *p=&y;
   cout<<"Please input two numbers:"<<endl;
   cin>>x;
   cin>>y;
   CompareIntNumberAndIntPointNumber1(x,p);
}
//编写一个函数，令其交换两个int指针
void ToChangePoint(int *x,int *y){
   int m=1;
   int n=2;
   x=&m;
   y=&n;
   cout<<"the address of point is :"<<x<<" "<<y<<endl;
   int *temp;
   temp=x;
   x=y;
   y=temp;   
   cout<<"the address changed of point is:"<<x<<" "<<y<<endl;
}
//函数递归调用实现阶乘
int DeGuiDiaoYong(int val){
   if (val>1)
   {
      int sum=1;
      sum=DeGuiDiaoYong(val-1)*val;
      return sum;
   }
   else
   {
      return 1;
   }
}
void DeGuiDiaoYongShow(){
   cout<<"Please input a numbers:"<<endl;
   int x;
   cin>>x;
   if (x>1)
   {
      cout<<DeGuiDiaoYong(x)<<endl;
   }
   else{
      cout<<"the number is valid"<<endl;
   }
}
void TheValuEqualXiaBiao(){
   int arr[10];
   for (int i = 0; i < 10; i++)
   {
      arr[i]=i;
   }
   for (int i = 0; i < 10; i++)
   {
      cout<<arr[i]<<" ";
   }
   cout<<endl;//换行
}
int TheValuEqualXiaBiao(int Arsize){
   int arr[Arsize];
   for (int i = 0; i < Arsize; i++)
   {
      arr[i]=i;
   }
   for (int *p=arr; *p!=arr[Arsize] ; p++)
   {
      cout<<*p<<" ";
   }
  
   return 0;
}

//利用begin和and函数做数组的遍历
void BeginAndToShowArr(){
   int arr[]={0,1,2,3,4,5};
   int *beg=begin(arr);
   int *last=end(arr);
   while (beg!=last)
   {
      cout<<*beg<<" ";
      beg++;
   } 
}
//求解数组的最大连续和
void MaxSumOfArr(){
   int total=0;
  int A[]={1,2,3,-1,5,6,7,8,9,1,-2};
  int best=A[1];
  int n=11;
  int sum1=0;
  for (int i = 0; i < n; i++)
  {
   for (int j = i; j < n; j++)
   {
     int sum=0;
     for (int k = i; k <=j ; k++)
     {
      sum=sum+A[k];
      total++;
     }
       if (sum>best)
      {
         best=sum;
      }
     
   }
   
  }
  cout<<best<<" "<<total<<endl;

  for (int i = 0; i < n; i++)
  {
   sum1=sum1+A[i];
  }
  cout<<sum1<<endl; 
}

//const的使用
void ConstUse(int& ref){
    ref=30;
}
void ConstShow(){
   int a=10;
   cout<<a<<endl;
   int &ref=a;
   ref=20;
   cout<<a<<endl;
   ConstUse(a);
   cout<<a<<endl;
}
//从任意输入的三个整数。找出按大小顺序排序处于中间位置的数。
void MidShow(int x,int y,int z){
   if ((y-x)*(y-z)<=0)
   {
      cout<<"the mid numbers is :"<<y<<endl;
   }
   else if ((x-y)*(x-z)<=0)
   {
      cout<<"the mid numbers is :"<<x<<endl;
   }
   else
   {
      cout<<"the mid numbers is :"<<z<<endl;
   }
     
}
//4.1.2逆序输出数组元素

void NiXuOutPutArr(int size){
   double arr[size]={};
   arr[0]=2;
   for (int i = 0; i < size; i++)
   {
      arr[i+1]=2*arr[i];
   }
   for (int i = 0; i < size; i++)
   {
      cout<<arr[size-1-i]<<" ";
      
   }
    for (int i = 0; i < size; i++)
   {
      cout<<arr[i]<<" ";
      
   }
}
//整数数组，n个元素，用键盘输入，将数组第一个元素移到数组末位，其余元素
//依次前移一个位置输出
void MoveArr(int size){
   int arr[size];
   for (int i = 0; i < size; i++)
   {
      arr[i]=i;
   }
   for (int i = 0; i < size; i++)
   {
      cout<<arr[i]<<" ";
   }
   for (int i = 0; i < size-1; i++)
   {
      int temp;
      temp=arr[i];
      arr[i]=arr[i+1];
      arr[i+1]=temp;
   }
   for (int i = 0; i < size; i++)
   {
      cout<<arr[i]<<" ";
   }
}
//分糖
void SharingCandy(int size){
   int arr[size]={};
   for (int i = 0; i < size; i++)
   {
      arr[i]=2*i+3;
   }
   for (int i = 0; i < size; i++)
   {
      cout<<i+1<<":"<<arr[i]<<" ";
   }
   for (int i = 0; i < size-1; i++)
   {
      arr[i+1]=arr[i]/2+arr[i+1];
      arr[i]=arr[i]/2;
   }
   cout<<"-----------"<<endl;
   for (int i = 0; i < size; i++)
   {
      cout<<i+1<<":"<<arr[i]<<" ";
   }
}
//灯光控制，二的倍数的灯变状态，三的倍数的灯改变状态
void LightChange(){
   int arrsize=1000;
   int arr[arrsize]={};
   //初始时所有灯熄灭
   for (int i = 0; i < arrsize; i++)
   {
      arr[i]=0;
   }
   for (int i = 0; i < arrsize; i++)
   {
      cout<<i+1<<":"<<arr[i]<<" ";
   }
   for (int i = 0; i < arrsize-1; i++)
   {
      if ((i+1)%2==0)
      {
         arr[i+1]=1;
         if ((i+1)%3==0)
         {
            arr[i+1]=abs(arr[i+1]-1);
         }  
         else
         {
            arr[i]=0;
         } 
      }
   }
   cout<<"-----------------"<<endl;
   for (int i = 0; i < arrsize; i++)
   {
       cout<<i+1<<":"<<arr[i]<<" ";
   }
   
}
//二维数组赋值输出
void TwoArrShow(int h,int l){
   
   int arr[h][l]={0};

    for (int i = 0; i < h; i++)
   {
      for (int j = 0; j < l; j++)
      {
        arr[i][j]=0; 
      }
   }

   for (int i = 0; i < h; i++)
   {
      for (int j = 0; j < l; j++)
      {
         cout<<arr[i][j]<<" ";
         
      }  
      cout<<"\n";
   }

   for (int i = 0; i < h; i++)
   {
      for (int j = 0; j < l; j++)
      {
        arr[i][j]=1; 
      }
   }
   for (int i = 0; i < h; i++)
   {
      for (int j = 0; j < l; j++)
      {
         cout<<arr[i][j]<<" ";
      }
      cout<<"\n";  
   }
}
//矩阵加法计算

   






//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//vim配置文件操作步骤
//1.在根目录/下进行操作  cd /etc/vim
//2.然后进行操作 sudo vim vimrc
//3.进行写操作
//4.将下列命令复制到命令末尾
/*
set nu                      " 显示行号
set tabstop=4               " 设置软制表符宽度为4
set softtabstop=4           " 设置软制表符宽度为4
set shiftwidth=4            " 设置缩进的空格数为4
set autoindent              " 设置自动缩进：即每行的缩进值与上一行相等
set cindent                 " 使用 C/C++ 语言的自动缩进方式
set cursorline              " 突出显示当前行
set expandtab               " 空格代替制表符
set showmatch               " 光标遇到圆括号、方括号、大括号时，自动高亮对应的另一个圆括号、方括号和大括号
set ruler                   " 在状态栏显示光标的当前位置（位于哪一行哪一列）
 
set guifont=Consolas:h15    " 设置字体和字体大小
 
set nobackup                " 取消备份文件
                            " 默认情况下，文件保存时，会额外创建一个备份文件，它的文件名是在原文件名的末尾，再添加一个波浪号~
setlocal noswapfile         " 不创建交换文件。交换文件主要用于系统崩溃时恢复文件，文件名的开头是.、结尾是.swp
set noundofile              " 取消生成un文件
 
set hlsearch                " 设置高亮显示搜索字符串
set showmode                " 在底部显示，当前处于命令模式还是插入模式
set showcmd                 " 命令模式下，在底部显示，当前键入的指令。比如输入快捷键将在底部显示具体命令
set t_Co=256                " 启用256色
set noerrorbells            " 出错时不要发出响声
" 高亮显示
syntax on
syntax enable 
*/
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//vim编辑器的操作键位
/*
vim模式
普通模式 NORMAL
插入模式 INSERT
可视模式 VISUAL
命令模式：
普通模式：
   在普通模式下，选中要复制的代码，然后按住ctrl，然后就能拖动复制选中代码
转换
普通模式 NORMAL  enter i (在光标前插入)> 插入模式 INSERT
普通模式 NORMAL  enter I (在行首插入)> 插入模式 INSERT
普通模式 NORMAL  enter a (在光标后插入)> 插入模式 INSERT
普通模式 NORMAL  enter A (在行尾插入)> 插入模式 INSERT
普通模式 NORMAL  enter o (在下一行插入)> 插入模式 INSERT
普通模式 NORMAL  enter O (在上一行插入)> 插入模式 INSERT
插入模式 INSET  enter ESC/jj/CapLocK > 普通模式NORMAL
普通模式 NORMAL enter ESC/v > 可视模式 VISUAL
普通模式 NORMAL enter : > 命令模式：
命令模式： enter ESC > 普通模式 NORMAL enter 
光标的移动
普通模式
k上，h左，j下，l右
w 跳到一个单词开头
b 跳到本单词或上一单词开头
e 跳到本单词或下一个单词结尾
ge 跳到上一个单词结尾
f{char} 光标跳到下一个{char}所在位置
F{char} 反向移动到上一个{char}所在位置
t{char} 光标跳到下一个{char}的前一个字符位置
T{cahr} 光标反向移动到上个{char}的后一个字符的位置
；重复上次的字符查找操作
，反向查找上次的查找命令
前面是操作符告诉vim要做什么
而动作motion是告诉vim我们要怎么做
操作符
d(delete)删除
c(change)修改（删除并进行插入模式）
y(yank)复制
p(paste)粘贴
d(撤销)
v(visual)选中并进入VISUAL模式
*/

