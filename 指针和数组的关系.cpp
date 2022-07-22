/*
 * @Author: ‘15071832337’ ‘1418868984@qq.com’
 * @Date: 2022-07-22 11:39:19
 * @LastEditors: ‘15071832337’ ‘1418868984@qq.com’
 * @LastEditTime: 2022-07-22 19:44:32
 * @FilePath: \vscode\1.cpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include <iostream>
 void a1();
 void a2();
 void a3();

using namespace std;
const int MAX = 4;

//定义获取数组长度 sizeof会报错

 
int main ()
{
   a3();
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



