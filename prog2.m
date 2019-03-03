function out=prog2(af)


af='/Users/Seven/Desktop/pic/';
addpath /Users/Seven/Desktop/siftDemoV4
S = dir(strcat(af,'*.JPG')); % pattern to match filenames.
disp(length(S));
r=zeros(3,3,length(S))
r(:,:,1)=[1 0 0;0 1 0;0 0 1]
yd=[0 0]
F = strcat(af,S(round(length(S)/2)).name);
    I1 = imread(F);
    [rr,cc,th]=size(I1);
    re=I1;
bound=[0 cc 0 rr]
disp(length(S));
for z =1:-2:-1
    yy=0;
    if z==-1 &length(S)/2==round(length(S)/2)
        yy=-1;
    end
for k = 1:floor(length(S)/2)+yy                            %
    F = strcat(af,S(round(length(S)/2)+z*(k-1)).name);
    I1 = imread(F);
    F = strcat(af,S(round(length(S)/2)+z*k).name);
    I2 =imread(F);
I=rgb2gray(I1);
II=rgb2gray(I2);
[rr,cc]=size(II);
%[num, matchidx, loc1, loc2] = match (strcat(af,S(k).name), strcat(af,S(k+1).name));
%img1_pts = loc1(find(matchidx > 0), 1:2);
%img2_pts = loc2(matchidx(find(matchidx > 0)), 1:2);
points1 = detectSURFFeatures(I);
points2 = detectSURFFeatures(II);
[f1,vpts1] = extractFeatures(I,points1);
[f2,vpts2] = extractFeatures(II,points2);
indexPairs = matchFeatures(f1,f2);
img1_pts = vpts1(indexPairs(:,1)).Location;
img2_pts = vpts2(indexPairs(:,2)).Location;
%figure; showMatchedFeatures(I,II,img1_pts,img2_pts);
imgp=[img1_pts img2_pts];
%disp(imgp.Location);
[num,num2]=size(imgp);
%disp(num);disp(num2);
v=0.0
while v<0.78
    rl=zeros(8,9) 
    rll=[];
    for i=0:3
        t=uint8(rand*(num-i));
        rll=[rll;imgp(t,:)]
        rl(2*i+1,:)=[-1.*imgp(t,3:4) -1 0 0 0 imgp(t,1).*imgp(t,3:4) imgp(t,1)];
        if i~=4
        rl(2*i+2,:)=[0 0 0 -1.*imgp(t,3:4) -1 imgp(t,2).*imgp(t,3:4) imgp(t,2)];
        end
        imgp(t,:)=[];
    end
    sb=imgp;
    imgp=[sb;rll];
    [a,b,c]=svd(rl,0);
    h=c(:,end);
    ht=h';
    hh=[ht(1:3);ht(4:6);ht(7:9)];
    %disp(h);
    e=0;
    e1=0;
    for i=1:num
        xx=hh*([img2_pts(i,:),1]')
        e1=abs(xx(1,1)/xx(3,1)-img1_pts(i,1))+abs(xx(2,1)/xx(3,1)-img1_pts(i,2));
        if e1>9
            e=e+1;
        end
    end
    v=1-e/num
end
%disp(r(:,:,k))
r(:,:,k+1)=r(:,:,k)*hh;
%disp(r(:,:,k+1));
c4=[0 0 rr rr;0 cc 0 cc;1 1 1 1]
rc=r(:,:,k+1)*c4;
for i=1:4
    if rc(1,i)/rc(3,i)<bound(1)
        bound(1)=round(rc(1,i)/rc(3,i));
    end
    if rc(1,i)/rc(3,i)>bound(2)
        bound(2)=round(rc(1,i)/rc(3,i));
    end
    if rc(2,i)/rc(3,i)<bound(3)
        bound(3)=round(rc(2,i)/rc(3,i));
    end
    if rc(2,i)/rc(3,i)>bound(4)
        bound(4)=round(rc(2,i)/rc(3,i));
    end
end

               % figure
              %imshow(I1)
%re=zeros(bound(4)-bound(3),bound(2)-bound(1),3);
[r1,c1,sbl]=size(re);
     tre=re;
     re=[zeros(r1,-bound(1)-yd(1),3) tre zeros(r1,bound(2)-c1+yd(1),3)]
     tre=re;
     [r1,c1,sbl]=size(re);
     %disp(size(re))
     re=[zeros(-bound(3)-yd(2),c1,3);tre;zeros(bound(4)-r1+yd(2),c1,3)]
     yd=[-bound(1) -bound(3)];
     %disp(size(re))
     %disp(yd)
for i=1:rr
    for j=1:cc
        nc=r(:,:,k+1)*[j;i;1];
        nc=nc./nc(3,1);
        nc1=r(:,:,k+1)*[j+1;i+1;1];
        nc1=nc1./nc1(3,1);
        iij=1;jjj=1;
        
                %disp(nc);disp(nc1);
                %disp(bound);
        if nc(1)>nc1(1)
            iij=-1;
        end
        if nc(2)>nc1(2)
            jjj=-1;
        end
        for ii=yd(1)+round(nc(1)):iij:yd(1)+round(nc1(1))-iij
            for jj=yd(2)+round(nc(2)):jjj:yd(2)+round(nc1(2))-jjj
                %disp(ii);disp(jj);
                if jj>0&ii>0
                re(jj,ii,:)=I2(i,j,:);
                end
            end
        end
    end
end
%hold on; 
figure;
imshow(re);
%plot(xs, ys, 'b-');
%drawnow;
                
end
end
     re=zeros(bound(4)-bound(3),bound(2)-bound(1),3);
     for i=rr
                
    %imshow(I)
    %S(k).data = I; % optional, save data.
end
%pa2("download/liptracking2/",1,4,2);
