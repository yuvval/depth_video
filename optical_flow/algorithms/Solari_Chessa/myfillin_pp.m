function O=myfillin_pp(Oin)
%2nd formula@p

filt = fspecial('gaussian', [15 15], 15/6);
filt=filt./sum(sum(fspecial('gaussian', [15 15], 15/6)));

tmp=Oin(:,:,1); tmp2=Oin(:,:,1);
tmp2(isnan(tmp2))=0;
masknonan=~isnan(tmp);
trueborder = bwmorph(masknonan,'remove');
tmp((trueborder)==0)=0;
num=conv2b(tmp, filt, 3);
den=conv2b(trueborder, filt, 3);
den(den==0)=1;
O1filled=(num./den).*(~masknonan) + tmp2.*masknonan;
clear tmp tmp2 masknonan trueborder num den

tmp=Oin(:,:,2); tmp2=Oin(:,:,2);
tmp2(isnan(tmp2))=0;
masknonan=~isnan(tmp);
trueborder = bwmorph(masknonan,'remove');
tmp((trueborder)==0)=0;
num=conv2b(tmp, filt, 3);
den=conv2b(trueborder, filt, 3);
den(den==0)=1;
O2filled=(num./den).*(~masknonan) + tmp2.*masknonan;
clear tmp tmp2 masknonan trueborder num den


O(:,:,1)=O1filled;
O(:,:,2)=O2filled;





