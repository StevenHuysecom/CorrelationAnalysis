function [row,col,e,MagX,MagY] = phasor(array)

fftROI = fft(array);
MagX = abs(fftROI(1,1));

end
