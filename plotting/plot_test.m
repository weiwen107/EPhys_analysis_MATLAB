ci = 1;

for ti = 34:48
    figure;
    plot(1:numel(nonzeros(width{1,ci}(ti,:))),nonzeros(width{1,ci}(ti,:)),'o')
end
