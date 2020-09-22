function [note] = notecreate(frq_no, stop)
note = sin(2*pi* [1:stop]/8192 * (440*2.^((frq_no-1)/12)));
end 