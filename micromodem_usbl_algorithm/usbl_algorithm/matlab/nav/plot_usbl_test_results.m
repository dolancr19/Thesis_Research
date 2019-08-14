
plot_timing = 0;
plot_std = 0;
plot_noise = 0;
plot_time = 1;

if plot_timing
  load usbl_timing.txt
  fs = 1/8e-6;
  ts = [1:length(usbl_timing)]./fs;
  figure(1);clf
  plot(usbl_timing(:,1), usbl_timing(:, 2), 'g')
  hold on
  plot(usbl_timing(:,1), usbl_timing(:, [6 5]), 'linewidth', 1.5)
  hold off
  legend('10 ms probe', 'USBL alg. time', 'alg+printing time')
  xlabel('seconds')
  ylabel('volts')
  title('Timing for USBL with 10 ms probe, 115k baud')
end

if plot_std
  nreps = 100;
  azsent =-90:10:90;
  elsent =-45;
  filename = 'um2ex100reps.txt';
  [dnum,az,el,owtt, irphase,t_ms] = read_usbl_csv(filename);
  az_reps = reshape(az, nreps, length(azsent)*length(elsent));
  el_reps = reshape(el, nreps, length(azsent)*length(elsent));
  az_er = az_reps - ones(nreps,1)*azsent;
  figure(2);clf;subplot(2,1,1);
  azstd = std(az_er);
  plot(azsent, azstd);grid on
  xlabel('azimuth, deg')
  ylabel('deg')
  title('std. deviation of error in azimuth at el=-45^0, alg. exhaustive Search')
  subplot(2,1,2);
  az_er_mean = mean(abs(az_er));
  plot(azsent, az_er_mean);grid on
  xlabel('azimuth, deg')
  ylabel('deg')
  title('Mean error at el=-45^0, alg. exhaustive search')
  orient tall
  print('-dpdf', 'um2ex100reps.pdf')
end

if plot_time
  nreps = 1;
  azsent =-90:10:90;
  azsent =-90:10:90;
  filename = 'um2gmtiming.txt';
  [dnum,az,el,owtt, irphase, t_ms] = read_usbl_csv(filename);
  figure(1);clf;subplot(2,1,1);
  hist(t_ms);
  xlabel('ms')
  tstr = sprintf('Golden mean search time, mean %dms, std %d', round(mean(t_ms)), round(std(t_ms)));
  title(tstr)
  subplot(2,1,2);
  filename = 'um2extiming.txt';
  [dnum,az,el,owtt, irphase, t_ms] = read_usbl_csv(filename);
  subplot(2,1,2);
  hist(t_ms);
  xlabel('ms')
  tstr = sprintf('Exhaustive search time, mean %dms, std %d', round(mean(t_ms)), round(std(t_ms)));
  title(tstr)
 orient tall
  print('-dpdf', 'um2timing.pdf')
end

