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
colorscheme molokai         " 设置主题为molokai
 
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



