function draw_lines(lines1, lines2, params);
% draws two groups of line segments

H = params.H;
W = params.W;
hfig = figure();
hold on
for i=1:size(lines1,1)
   l = lines1(i,:);
    plot([l(1) l(3)], H - [l(2) l(4)]);
end

for i=1:size(lines2,1)
   l = lines2(i,:);
    plot([l(1) l(3)], H - [l(2) l(4)],'r');
end

axis([0 W 0 H])
axis off
axis equal
print(hfig,'-depsc', '-r300', sprintf('%s/lines', params.folder_out));