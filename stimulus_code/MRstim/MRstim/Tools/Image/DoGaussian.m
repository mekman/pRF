function [final,varargout]=DoGaussian(im, scale)
% DoGaussian - gaussian blurring of image
%  [final,varargout]=DoGaussian(im, scale)
%   scale is standard dev of Gaussian

% from Steven Dakin

NoFilters=length(scale);
[m n,nchannels]=size(im);
if (nchannels>1)
    final=zeros(m,n,nchannels,NoFilters);
else
    final=zeros(m,n,NoFilters);
end
for i=1:NoFilters
    if (scale(i)>0)
        log1(:,:,i) = fft2(fspecial('gaussian',[m n],scale(i)));
    end
end
for i=1:nchannels
	if (isreal(im(:,:,i)))
			ffttmp(:,:,i)=fft2(double(im(:,:,i)));
		else
			ffttmp(:,:,i)=double(im(:,:,i));
        end
 end
 for i=1:nchannels
    for j=1:NoFilters
        if (scale(j)>0)
			ThisRes=real(fftshift(real(ifft2(ffttmp(:,:,i).*log1(:,:,j)))));
        else
            ThisRes=im;
        end            
        if (nchannels>1)
               final(:,:,i,j)=ThisRes; 
        else
               final(:,:,j)=ThisRes; 
        end

	end;		
end
if nargout>1
    for i=1:nchannels
    for j=1:NoFilters
                if (nchannels>1)
               a=final(:,:,i,j);
               energy1(i,j)=std(a(:));
        else
               a=final(:,:,j); 
               energy1(j)=std(a(:));
        end
    end
    end 
   varargout={energy1};
end
