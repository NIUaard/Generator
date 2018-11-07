function Cropped_Image=Crop_Image (A);

fprintf('please choose the top left and bottom right for the part of the image\n')
imagesc(A)
%axis xy
MM= ginput(2);
MinX = round(MM(1,1));
MaxX = round(MM(2,1));
MinY = round(MM(1,2));
MaxY = round(MM(2,2));

MM


Cropped_Image = A(MinY:MaxY,MinX:MaxX);


