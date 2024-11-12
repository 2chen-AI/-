Consider the third order plant

$$y=G(x)u$$

where

$$G(s)=\dfrac{b_2s^2+b_1s+b_0}{s^3+a_2s^2+a_1s+a_0}$$

(a) For the case where all parameters are unknown, obtain a static linear parametric model (SPM) of the plant. Assume that $u$ and $y$ are available for measurement, but their derivatives are not.

(b) For the case where it is known that $$a_0=b_2=2$$, obtain a parametric model for the plant in terms of $$\theta^*=[b_1,b_0,a_2,a_1]^T$$.

(c\) For the case where it is known that $$b_0=b_1=0, b_2=3$$ obtain a parametric model in terms of $$\theta^*=[a_2,a_1,a_0]^T$$. Design and simulate in Matlab/Simulink a gradient based algorithm to generate on line esrimates of $\theta^*$. For the actual system you can take $a_0=b_2=2,b_0=b_1=0,b_2=3,a_1=a_2=2$. First use $u(t)=1$, then use $u(t)=\sin(0.2t)+\sin(t)$. Compare the different results.
Solution:
(a) since $u$ and $y$ are available for measurement, but their derivatives are not,
both sides time $s^3 + a_2s^2+a_1s+a_0$, then
$$
(s^3+a_2s^2+a_1s+a_0)y=(b_2s^2+b_1s+b_0)u
$$solution 1
$$y'''+a_2y''+a_1y'+a_0y=b_2u''+b_1u'+b_0u$$
then  we can obtain a static linear parametric model (SPM) of the plant
$$y'''={\theta^*}^T\phi$$
where$${\theta^*}=\begin{bmatrix}a_2\\ a_1\\ a_0\\ b_2\\ b_1\\ b_0
\end{bmatrix},\phi=\begin{bmatrix}-y''\\ -y'\\ -y\\ u''\\ u'\\ u
\end{bmatrix}$$
emit the term containing derivatives, then we get $a_0y=b_0u$
So the static linear parametric mode of the plant is 
$$
y=\dfrac{b_0}{a_0}u
$$












(b) when $a_0=b_0=2$
$G(s)=\dfrac{2s^2+b_1s+b_0}{s^3+a_2s^2+a_1s+2}$
$$y=\dfrac{2s^2+b_1s+b_0}{s^3+a_2s^2+a_1s+2}u$$
$$(s^3+a_2s^2+a_1s+2)y=(2s^2+b_1s+b_0)u$$
$$\begin{equation*}\begin{split}&\underbrace{\dfrac{s^3}{s^3+\lambda_2s^2+\lambda_1s+\lambda_0}y}_{z}+a_2\dfrac{s^2}{s^3+\lambda_2s^2+\lambda_1s+\lambda_0}y+a_1\dfrac{s}{s^3+\lambda_2s^2+\lambda_1s+\lambda_0}y+2\dfrac{1}{s^3+\lambda_2s^2+\lambda_1s+\lambda_0}y\\-&2\dfrac{s^2}{s^3+\lambda_2s^2+\lambda_1s+\lambda_0}u-b_1\dfrac{s}{s^3+\lambda_2s^2+\lambda_1s+\lambda_0}u-b_0\dfrac{1}{s^3+\lambda_2s^2+\lambda_1s+\lambda_0}u=0\end{split}\end{equation*}$$
then we can obtain a static linear parametric model (SPM) of the plant
$$z={\theta^*}^T\phi$$
where
$$\theta^*=\begin{bmatrix}a_2\\ a_1\\ 2\\ 2\\ b_1\\ b_0
\end{bmatrix},\phi=\begin{bmatrix}-\frac{s^2}{s^3+\lambda_2s^2+\lambda_1s+\lambda_0}y\\ -\frac{s}{s^3+\lambda_2s^2+\lambda_1s+\lambda_0}y\\ -\frac{1}{s^3+\lambda_2s^2+\lambda_1s+\lambda_0}y\\ \frac{s^2}{s^3+\lambda_2s^2+\lambda_1s+\lambda_0}u\\ \frac{s}{s^3+\lambda_2s^2+\lambda_1s+\lambda_0}u\\ \frac{1}{s^3+\lambda_2s^2+\lambda_1s+\lambda_0}u
\end{bmatrix}$$
emit the term containing derivatives, then we get
$$
2y=b_0u
$$considering $b_1,a_2,a_1$,we get
$$
y=\dfrac{b_0}{2}u+\dfrac{b_1}{2}\cdot\dfrac{du}{dt}-\dfrac{a_2}{2}\cdot\dfrac{dy}{dt}-\dfrac{a_1}{2}\cdot \dfrac{dy}{dt}
$$so
$$
y=\dfrac{1}{2}{\theta^*}^T[\dfrac{du}{dt},u,-\dfrac{dy}{dt},-\dfrac{dy}{dt}]
$$


















(c\) when $b_0=b_1=0,b_2=3$
$$G(s)=\dfrac{3s^2}{s^3+a_2s^2+a_1s+a_0}$$
then
$$y=\dfrac{3s^2}{s^3+a_2s^2+a_1s+a_0}u$$
$$(s^3+a_2s^2+a_1s+a_0)y=3s^2u$$
$$\begin{equation*}
\begin{split}&
\underbrace{\dfrac{s^3}{s^3+\lambda_2s^2+\lambda_1s+\lambda_0}y}_z+a_2\dfrac{s^2}{s^3+\lambda_2s^2+\lambda_1s+\lambda_0}y\\+&a_1\dfrac{s}{s^3+\lambda_2s^2+\lambda_1s+\lambda_0}y+a_0\dfrac{1}{s^3+\lambda_2s^2+\lambda_1s+\lambda_0}y-3\dfrac{s^2}{s^3+\lambda_2s^2+\lambda_1s+\lambda_0}u=0
\end{split}
\end{equation*}$$

then we can obtain a parametric model
$$z={\theta^*}^T\phi$$
where
$$\theta^*=\begin{bmatrix}a_2\\ a_1\\ a_0\\ 3
\end{bmatrix},\phi=\begin{bmatrix}-\frac{s^2}{s^3+\lambda_2s^2+\lambda_1s+\lambda_0}y\\ -\frac{s}{s^3+\lambda_2s^2+\lambda_1s+\lambda_0}y\\ -\frac{1}{s^3+\lambda_2s^2+\lambda_1s+\lambda_0}y\\ -\frac{s^2}{s^3+\lambda_2s^2+\lambda_1s+\lambda_0}u
\end{bmatrix}$$


Design and simulate in Matlab a gradient based algorithm to generate on line estimates of $\theta^*$. For the actual system we take $a_0 =2, b_0 = b_1 =
 0, b_2 = 3,a_1 = a_2 = 2$.  First use $u(t) = 1$, then use $u(t) = \sin (0.2t) + \sin (t)$. Compare the different results. 
```matlab
// An highlighted block
clear,clc
% 定义系统真实参数
a0_true = 2; a1_true = 2; a2_true = 2;

% 初始化参数估计值
theta_hat = [1; 1; 1];

alpha = 0.01;% 定义学习率
epsilon = 1e-6; % 小扰动值

% 定义模拟时间步长和总时间
dt = 0.1; T = 100; t = 0:dt:T;

% 定义输入信号
% u = ones(size(t));%u(t)=1
u = sin(0.2*t) + sin(t);%
u = u(:); % 将输入信号转换为列向量

% 初始化输出估计值和误差
y_hat = zeros(size(t));
e = zeros(size(t));

% 根据实际参数计算整个时间序列的估计输出
G_true = tf([3 0 0], [1 2 2 2]);
y_true = lsim(G_true, u, t); % 一次性计算整个时间序列的估计输出

for k = 1:1
    % 根据当前估计参数计算整个时间序列的估计输出
    G_hat = tf([3 0 0], [1 theta_hat(1) theta_hat(2) theta_hat(3)]);
    y_hat = lsim(G_hat, u, t); % 一次性计算整个时间序列的估计输出
    
    % 计算整个时间序列的误差
    e = y_true - y_hat; % 这里假设实际输出为y_true，根据系统传递函数关系
    
    % 计算梯度（使用有限差分法）
    J_theta2_original = (e.^2)/2; % 损失函数，误差的平方
    grad_J = zeros(size(theta_hat));
    for i = 1:length(theta_hat)
        theta_hat_perturbed = theta_hat;
        theta_hat_perturbed(i) = theta_hat_perturbed(i) + epsilon;
        G_hat_perturbed = tf([3 0 0], [1 theta_hat_perturbed(1) theta_hat_perturbed(2) theta_hat_perturbed(3)]);
        y_hat_perturbed = lsim(G_hat_perturbed, u, t);
        e_perturbed = y_true - y_hat_perturbed;
        J_theta2_perturbed = (e_perturbed.^2)/2;
        grad_J(i) = mean((J_theta2_perturbed - J_theta2_original) / epsilon); % 求平均后赋值
    end
    
    % 更新参数估计值
    theta_hat = theta_hat - alpha * grad_J;
end
% 比较结果（可以根据需要进一步分析和可视化结果）
% disp(['使用u(t)=1时最终参数估计值: ', num2str(theta_hat(1)), ' ', num2str(theta_hat(2)), ' ', num2str(theta_hat(3))]);
disp(['使用u(t)=sin(0.2t)+sin(t)时最终参数估计值: ', num2str(theta_hat(1)), ' ', num2str(theta_hat(2)), ' ', num2str(theta_hat(3))]);
```
