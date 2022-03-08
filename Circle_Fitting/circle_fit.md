# 利用散点估计外接圆


一般而言，对散点进行圆弧拟合采用最小二乘方法，而最小二乘不可避免的需要求解线性方程组，且矩阵的维度随着散点的数量增加而增加，这使得在嵌入式设备上用于估计道路的一定范围的曲率较为耗费算力。

那么肯定会有人说，如果我基于散点得到道线的解析描述，这样不就可以得到求出曲率，不这样做的原因有以下两点：

   -  这样求出曲率是逐点的曲率，在实际问题中一定范围的平均曲率则更为实用 
   -  求出的道线的多项式无论是基于 Hough 变换还是最小二乘本身是有一定的误差，而且由于求曲率需要涉及到高阶导数，会使得拟合的误差对于结果影响很大 

所以下面介绍一种在一定前提下对于最小二乘法的简化方法

# 算法原理

圆方程可以写成：

$$
(x-x_c)^2 + (y - y_c)^2 = R^2
$$
那么散点拟合圆的距离平方和就要求最小：

$$
min \sum_{i}^n (\sqrt{(x_i-x_c)^2 + (y_i - y_c)^2} - R)^2
$$
注意，上式难以直接求解，需要形成矩阵方程来求解。

**为了简化，我们退而求其次，求解下式**
$$
min \sum_{i}^n ((x_i-x_c)^2 + (y_i - y_c)^2 - R^2)^2
$$
本质上是将原式展开并采用一个很强的假设：

<img src="https://latex.codecogs.com/gif.latex?\sqrt{(x_i&space;-x_c&space;)^2&space;+(y_i&space;-y_c&space;)^2&space;}-R=\delta&space;R\to&space;0~~\forall&space;i\le&space;n"/>

显然，这个假设当采样点都几乎分布在圆弧的附近时才会成立，所以这个方法会非常容易受到离群点的影响。

# 推导

那么我们进一步可以定义一个辅助函数：

<img src="https://latex.codecogs.com/gif.latex?g(x_c&space;,y_c&space;,R)=(x_i&space;-x_c&space;)^2&space;+(y_i&space;-y_c&space;)^2&space;-R^2"/>

问题就转变为

<img src="https://latex.codecogs.com/gif.latex?min\sum_i^n&space;g(x_c&space;,y_c&space;,R)^2"/>

为无约束优化问题，其最值必要条件为各方向偏导数为0，分别有

<img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;f}{\partial&space;R}=-2R\times&space;\sum&space;\left({\left(x_i&space;-x_c&space;\right)}^2&space;+{\left(y_i&space;-y_c&space;\right)}^2&space;-R^2&space;\right)=-2R\times&space;\sum&space;g\left(x_i&space;,y_i&space;\right)=0"/>

显然圆的半径不为 0 ，则有 

<img src="https://latex.codecogs.com/gif.latex?\sum&space;g\left(x_i&space;,y_i&space;\right)=0"/>

进一步从对 <img src="https://latex.codecogs.com/gif.latex?\inline&space;x_c&space;~~y_c"/> 导数可以得到

<img src="https://latex.codecogs.com/gif.latex?\begin{array}{l}&space;\sum&space;x_i&space;g\left(x_i&space;,y_i&space;\right)=0\\&space;\sum&space;y_i&space;g\left(x_i&space;,y_i&space;\right)=0&space;\end{array}"/>

设置

<img src="https://latex.codecogs.com/gif.latex?\begin{array}{l}&space;u_i&space;=x_i&space;-x\\&space;u_C&space;=x_C&space;-\bar{x}&space;\\&space;v_i&space;=y_i&space;-\bar{y}&space;\\&space;v_C&space;=y_c&space;-\bar{y}&space;&space;\end{array}"/>

其中

<img src="https://latex.codecogs.com/gif.latex?\begin{array}{l}&space;\bar{x}&space;=\sum&space;x_i&space;/N\\&space;\bar{y}&space;=\sum&space;y_i&space;/N&space;\end{array}"/>

那么有

<img src="https://latex.codecogs.com/gif.latex?\begin{array}{l}&space;\sum&space;u_i&space;g\left(u_i&space;,v_i&space;\right)=0\\&space;\sum&space;v_i&space;g\left(u_i&space;,v_i&space;\right)=0&space;\end{array}"/>

展开有

<img src="https://latex.codecogs.com/gif.latex?\begin{array}{l}&space;\sum&space;\left(u_i^3&space;-2u_i^2&space;u_c&space;+u_i&space;u_c^2&space;+u_i&space;v_i^2&space;-2u_i&space;v_i&space;v_c&space;+u_i&space;v_c^2&space;-u_i&space;R^2&space;\right)=0\\&space;\sum&space;\left(u_i^2&space;v_i&space;-2u_i&space;v_i&space;u_c&space;+v_i&space;u_c^2&space;+v_i^3&space;-2v_i^2&space;v_c&space;+v_i&space;v_c^2&space;-v_i&space;R^2&space;\right)=0&space;\end{array}"/>

带入 <img src="https://latex.codecogs.com/gif.latex?\inline&space;\sum&space;u_i&space;=0~~\sum&space;v_i&space;=0"/>

<img src="https://latex.codecogs.com/gif.latex?\begin{array}{l}&space;\sum&space;\left(u_i^3&space;-2u_i^2&space;u_C&space;+u_i&space;v_i^2&space;-2u_i&space;v_i&space;v_C&space;\right)=0\\&space;\sum&space;\left(u_i^2&space;v_i&space;-2u_i&space;v_i&space;u_C&space;+v_i^3&space;-2v_i^2&space;v_C&space;\right)=0&space;\end{array}"/>

进一步定义

<img src="https://latex.codecogs.com/gif.latex?\begin{array}{l}&space;S_{uuu}&space;=\sum&space;u_i^3&space;\\&space;S_{vvv}&space;=\sum&space;v_i^3&space;\\&space;S_{uu}&space;=\sum&space;u_i^2&space;\\&space;S_{vv}&space;=\sum&space;v_i^2&space;\\&space;S_{uv}&space;=\sum&space;u_i&space;v_i&space;\\&space;S_{uuv}&space;=\sum&space;u_i^2&space;v_i&space;\\&space;S_{uvv}&space;=\sum&space;u_i&space;v_i^2&space;&space;\end{array}"/>

则上面的求和式可以写成

<img src="https://latex.codecogs.com/gif.latex?\begin{array}{l}&space;S_{uu}&space;u_C&space;+S_{uv}&space;v_C&space;=\frac{S_{uuu}&space;+S_{uvv}&space;}{2}\\&space;S_{uv}&space;u_C&space;+S_{vv}&space;v_C&space;=\frac{S_{uuv}&space;+S_{vvv}&space;}{2}&space;\end{array}"/>

则可以求出

<img src="https://latex.codecogs.com/gif.latex?\begin{array}{l}&space;u_C&space;=\frac{S_{uuv}&space;S_{uv}&space;-S_{uuu}&space;S_{vv}&space;-S_{uvv}&space;S_{vv}&space;+S_{uv}&space;S_{vvv}&space;}{2\left(S_{uv}^2&space;-S_{uu}&space;S_{vv}&space;\right)}\\&space;v_C&space;=\frac{-S_{uu}&space;S_{uuv}&space;+S_{uuu}&space;S_{uv}&space;+S_{uv}&space;S_{uvv}&space;-S_{uu}&space;S_{vvv}&space;}{2\left(S_{uv}^2&space;-S_{uu}&space;S_{vv}&space;\right)}\\&space;R^2&space;=\sum&space;\left({\left(x_i&space;-x_c&space;\right)}^2&space;+{\left(y_i&space;-y_c&space;\right)}^2&space;\right)&space;\end{array}"/>

# Matlab 实现

```matlab:Code
% Sample Generation
clear; 
R = 4 ;

x = -0.6:0.2:3;
y = sqrt(R^2 - x.^2);

[Flag,R,Xc,Yc] = circleFit(x,y);

plot(x,y);
hold on;
scatter(Xc,Yc);
```

```matlab:Code
function [Flag,R,Xc,Yc] = circleFit(X_in,Y_in)

    
    sum_x = 0;
    sum_y = 0;
    sum_x2 = 0;
    sum_y2 = 0;
    sum_x3 = 0;
    sum_y3 = 0;
    sum_xy = 0;
    sum_x1y2 = 0;
    sum_x2y1 = 0;
    
    
    
    N = length(X_in); 
    
    for i = 1 : N
        x = X_in(i);
        y = Y_in(i);
        x2 = x*x;
        y2 = y*y;
        sum_x = sum_x + x;
        sum_y = sum_y + y;
        sum_x2 = sum_x2 + x2;
        sum_y2 = sum_y2 + y2;
        sum_x3 = sum_x3 + x2*x;
        sum_y3 = sum_y3 + y2*y;
        sum_xy = sum_xy + x*y;
        sum_x1y2 = sum_x1y2 + x*y2;
        sum_x2y1 = sum_x2y1 + x2*y;
    end
    
    C = N * sum_x2 - sum_x * sum_x;
    D = N * sum_xy - sum_x * sum_y;
    E = N * sum_x3 + N * sum_x1y2 - (sum_x2 + sum_y2) * sum_x;
    G = N * sum_y2 - sum_y * sum_y;
    H = N * sum_x2y1 + N * sum_y3 - (sum_x2 + sum_y2) * sum_y;
    if abs(C * G - D * D) < 10^-6 || abs((D * D - G * C)) < 10^-6
        
        Xc = 0;
        Yc = 0;
        R = 0;
        Flag = uint8(0);
       
    else
        a = (H * D - E * G) / (C * G - D * D);
        b = (H * C - E * D) / (D * D - G * C);
        c = -(a * sum_x + b * sum_y + sum_x2 + sum_y2) / N;
    
        Xc = -0.5*c;
        Yc = -0.5*b;
        R = sqrt(a*a + b*b - 4*c)/2;
        Flag = uint8(1);
    end
    
end
```
