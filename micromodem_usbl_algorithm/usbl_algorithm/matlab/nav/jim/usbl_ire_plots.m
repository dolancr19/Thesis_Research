function h_vec = usbl_ire_plots(ire, peak_sample, title_str, xlim_vec, ylim_vec)

subplot(2,1,1);
hp1a=plot(real(ire));
hold on;
hp1b=plot(peak_sample, real(ire(peak_sample,:)), 'o');
hold off;
grid on;
hx1=xlabel('IRE sample number');
hy1=ylabel('Re(IRE)');
ht1=title(title_str);
xlim(xlim_vec);
ylim(ylim_vec);

subplot(2,1,2);
hp2a=plot(imag(ire));
hold on;
hp2b=plot(peak_sample, imag(ire(peak_sample,:)), 'o');
hold off;
grid on;
hx2=xlabel('IRE sample number');
hy2=ylabel('Im(IRE)');
ht2=title(title_str);
xlim(xlim_vec);
ylim(ylim_vec);

h_vec = [hp1a; hp1b; hx1; hy1; ht1; hp2a; hp2b; hx2; hy2; ht2];
return;
