t = (0:0.001:1)'; 
y = sin(2*pi*t);
ns = 0.04; %strength of noise
yn = y + ns*randn(size(t));

ynh = hilbert(yn); %hilbert transform of yn
instphase = atan(imag(ynh) ./ real(ynh)) %instantaneous phase of yn
plot(t, instphase) %plot the instantaneous phase