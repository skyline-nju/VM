# 用PYTHON封装VM的c++代码
先写好接口文件VMpy.i，然后在cmd执行
```
$ swig -python -c++ VMpy.i
```
得到VMpy_wrap.cxx和VMpy.py文件。用VS新建一个项目，windows桌面-->windows桌面向导-->动态链接库。添加相关的源文件、头文件以及接口文件。

**环境变量**

1. **PYTHON_INCLUDE** C:\Users\user\Miniconda3\include

2. **PYTHON_LIB** C:\Users\user\Miniconda3\libs\python36.lib

3. **NUMPY_INCLUDE** C:\Users\user\Miniconda3\Lib\site-packages\numpy\core\include\

**VS编译设置：**

1. 将PYTHON_INCLUDE和PYTHON_LIB分别添加到*项目属性*-->*VC++目录*下的*包含目录*和*库目录*。如果用到numpy, 需要在*包含目录*中添加NUMPY_INCLUDE.
2. 在*C/C++*-->*预处理器*-->*预处理器定义*中添加**SWIG_EXPORTS**。
3. 在*连接器*-->*输入*-->*附加依赖项*中添加*PYTHON_LIB*.
4. 在*常规*中将*目标文件名*改为*_$(ProjectName)*，*扩展名*改为*.pyd*。

生成解决方案，将这一文件拷贝到VMpy.py所在文件夹，然后就能在其他python程序中import VMpy.
