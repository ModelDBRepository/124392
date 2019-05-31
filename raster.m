function raster(fn)
% RASTER(FN, N) make raster plot of (N randomly chosen) cells in the
% AP file FN

if nargin==1, numproc = 1; end

fd = fopen(fn);
ap = textscan(fd, '%n%n%q');
fclose(fd);

if isempty(ap)
	fprintf(1,'Empty input\n');
	return
end
cells = unique(ap{1});
if nargin>=3
    indx = randperm(length(cells));
    indx = indx([1:N]);
else
    indx = 1:length(cells);
end

newplot
set(gcf,'Name','Spike Raster Plot', ...
    'NumberTitle','off', ...
    'toolbar', 'figure')

h = zeros(1,max(cells));
for i=1:length(indx)
    j = cells(indx(i));
    X = ap{2}(ap{1}==j);
    Y = j*ones(1,length(X));
	h(j+1) = line(X, Y, ...
		'marker', 'o', ...
		'linestyle', 'none', ...
		'markersize', 4, ...
		'markerfacecolor', 'black');
	xlabel('Time (s)')
	ylabel('cell number')
end
grid on

if length(ap)>=3 	  
	c = get(gcf, 'DefaultAxesColorOrder');
	maxc = size(c, 1);
	[cellnames cindx] = unique(ap{3});
	lstr = {};
	hf = [];
	for i=1:length(cindx)
		lstr{i} = cellnames{i};
		hf(i) = h(ap{1}(cindx(i))+1);
	end
	indxx = find(hf~=0);
	legend(hf(indxx), lstr{indxx});
	for i=1:length(indxx)
		indzz = strcmp(ap{3},lstr{indxx(i)})==1;
		lcells = unique(ap{1}(indzz));
		col = c(mod(i-1, maxc)+1, :);
		set(h(lcells+1), 'Color', col);
	end
end

uicontrol('style', 'push', ...
	'string', 'refresh', ...
	'FontSize',8, 'FontWeight', 'bold', ...
	'units', 'normal', ...
	'pos', [0.9 .01 0.08 .04], ...
	'backgroundcolor', 'black', ...
	'foregroundcolor', 'white', ...
	'call', ['raster(''', fn, ''',' num2str(numproc), ')']);

