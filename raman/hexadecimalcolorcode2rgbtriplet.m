function [rgb] = hexadecimalcolorcode2rgbtriplet(hex)

rgb = NaN(length(hex),3);
for j=1:length(hex)
    for i=1:3
        rgb(j,i) = hex2dec(extractBetween(hex(j),i*2,i*2+1))/255;
    end
end

end

