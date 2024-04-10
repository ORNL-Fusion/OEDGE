plot_triangles
hold on
for i=1:ntri
    text(knot(i).R,knot(i).Z,num2str(i));
%     pause
end