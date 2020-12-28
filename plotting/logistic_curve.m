x = -6:0.01:6;
y = 1./(1+exp(-x));
target_color = 1/255*[143,188,143];

figure;

hold on
rectangle('Position',[-5.9,0.45,11.9,0.1],...
    'Curvature',[0.1,0.1],...
    'EdgeColor',target_color,...
    'FaceColor',target_color);
 plot(x,y,...
    'Linewidth',3);
 yticks(0:0.2:1);
 xticks([]);
 
 ax = gca;
 ax.XAxis.LineWidth = 2;
 ax.YAxis.LineWidth = 2;