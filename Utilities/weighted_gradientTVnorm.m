% function [tvg] =weighted_gradientTVnorm2A(img,del)
function [tvg] =weighted_gradientTVnorm2A(img,weightx,weighty,weightz)

tvg=zeros(size(img));

Gx=diff(img,1,1);
Gy=diff(img,1,2);
Gz=diff(img,1,3);


Gx=cat(1,zeros(size(Gx(1,:,:))),Gx);
Gy=cat(2,zeros(size(Gy(:,1,:))),Gy);
Gz=cat(3,zeros(size(Gz(:,:,1))),Gz);

% for plane=2:size(img,3)
%     for row= 2:size(img,1)
%         for col=2:size(img,2)
%             
%             weightx(row,col,plane)=exp(-(Gx(row,col,plane)./del)^2);
%             weighty(row,col,plane)=exp(-(Gy(row,col,plane)./del)^2);
%             weightz(row,col,plane)=exp(-(Gz(row,col,plane)./del)^2);
% 
%         end
%     end
% end
nrm=weighted_safenorm(Gx,Gy,Gz,weightx,weighty,weightz);

tvg(1:end,1:end,1:end)= tvg(1:end,1:end,1:end)+((weightx(1:end,1:end,1:end).*Gx(1:end,1:end,1:end))+(weighty(1:end,1:end,1:end).*Gy(1:end,1:end,1:end))+(weightz(1:end,1:end,1:end).*Gz(1:end,1:end,1:end)))./nrm(1:end,1:end,1:end);
tvg(2:end-1,:,:)=tvg(2:end-1,:,:)-(weightx([2:end-1]+1,:,:).*Gx([2:end-1]+1,:,:))./nrm([2:end-1]+1,:,:);
tvg(:,2:end-1,:)=tvg(:,2:end-1,:)-(weighty(:,[2:end-1]+1,:).*Gy(:,[2:end-1]+1,:))./nrm(:,[2:end-1]+1,:);
tvg(:,:,2:end-1)=tvg(:,:,2:end-1)-(weightz(:,:,[2:end-1]+1).*Gz(:,:,[2:end-1]+1))./nrm(:,:,[2:end-1]+1);




% for plane=1:size(img,3)
%     for row= 1:size(img,1)
%         for col=1:size(img,2)
%           
%             tvg(row,col,plane)=tvg(row,col,plane)+((weightx(row,col,plane)*Gx(row,col,plane))+(weighty(row,col,plane)*Gy(row,col,plane))+(weightz(row,col,plane)*Gz(row,col,plane)))./nrm(row,col,plane);
%    
% %-----------------------------------------------------------------------------------------------            
%             if row <size(img,1)
%                 term1=(-1*weightx(row+1,col,plane)*Gx(row+1,col,plane))./...
%                 nrm(row+1,col,plane);
%                 tvg(row,col,plane)=tvg(row,col,plane)+term1;
%             end
% %-----------------------------------------------------------------------------------------------     
%                    
%             if col <size(img,2)
%                 term2=(-1*weighty(row,col+1,plane)*Gy(row,col+1,plane))./...
%                 nrm(row,col+1,plane);
%                 tvg(row,col,plane)=tvg(row,col,plane)+term2;
%             end
% 
% %-----------------------------------------------------------------------------------------------              
%             if plane <size(img,3)
%                 term3=(-1*weightz(row,col,plane+1)*Gz(row,col,plane+1))./...
%                 nrm(row,col,plane+1);
%                 tvg(row,col,plane)=tvg(row,col,plane)+term3;
%             end
%      
%     
%             
%         end
%     end
% end




end

function nrm=weighted_safenorm(Gx,Gy,Gz,weightx,weighty,weightz)

ksi=0.00000001;
nrm=sqrt(ksi+(weightx.*Gx.^2)+(weighty.*Gy.^2)+(weightz.*Gz.^2));
nrm(nrm==0)=1;

end

