% Illustrate step response of the second order system that was identified. 
% $A=-0.5787$ and  $B= 0.1572$. 
d = 14
A = -0.5787
B = 0.1572
sys = tf(B, [1 -A 0], 'inputdelay', d)
dsys = tf(B, [1 -A ], 'inputdelay', d)
figure
subplot(121)
[Y, T] = step(sys)
plot(T, Y, 'linewidth', 2, 'color', [0 0.4470 0.7410])
xlabel('Time [days]')
ylabel('Feedback measure z(t)')
xlim([0 50]); grid on;
subplot(122)
[Yd, Td] = step(dsys)
plot(Td, Yd, 'color', [0.8500 0.3250 0.0980], 'linewidth', 2)
xlabel('Time [days]')
ylabel('Derivative dz(t)/dt')
xlim([0 50]); grid on; 
set(gcf, 'color', [1 1 1])