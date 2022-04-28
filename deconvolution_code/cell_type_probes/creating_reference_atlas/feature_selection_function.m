%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Portions of this file were adapted from Joshua Moss and Tommy Kaplan from the paper                                %
% Comprehensive human cell-type methylation atlas reveals origins of circulating cell-free DNA in health and disease %
% doi: https://doi.org/10.1038/s41467-018-07466-6                                                                    %
% Any reuse of this code should include that attribution                                                             %               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% For Research Use only %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function feature_selection(ncpgs)
% read atlas
A = importdata('Atlas.X1.LUMP0.7.csv',',',1);
hdrs = A.textdata(1,2:end);
HD=hdrs; for i=1:length(HD), HD{i}(HD{i}=='_')=' '; end;
cg = A.textdata(2:end,1);
A = A.data;

[N,d] = size(A);
fprintf('Atlas loaded (1/2) [%dx%d]\n',N,d);

%	STD = importdata('STD.X1.LUMP0.7.csv',',',1); STD = STD.data;

% read CpG information
[~,chr,pos] = textread('CpGs.Illumina_450k.tab','%s%s%d%*[^\n]','headerlines',1,'delimiter','\t');
chrs = unique(chr); for i=1:length(chrs), c.(chrs{i}) = ismember(chr,chrs{i}); end; clear chrs;
chrs = unique(chr); [~,J]=ismember(chr,chrs); chr=J;

% npos - number of *on-array* neighboring CpGs
[~,npos] = textread('CpGs.Illumina_450k_blocks.tab','%s%d%*[^\n]','headerlines',1,'delimiter','\t');

fprintf('Data loaded\n');

%%%%%%%%%%%%%%%%%%%%%
% Feature selection
%%%%%%%%%%%%%%%%%%%%%

% ncpgs = 4000; 
flank=50; nECC=500; nECCmax=500; opts = struct('Accy',0);
fprintf('Searching for 2x%d CpGs, with %dbp flanks, with additional %d pairwise CpGs (%d per pair)\n', ...
	ncpgs, flank, nECC, nECCmax);

if 1,
    rng(6428049)
    I=find(var(A,[],2)>=0.001 & ~any(isnan(A),2));
    % I=find(var(A,[],2)>=0.001 & ~any(isnan(A),2) & max(STD,[],2)<=0.15);
    if ncpgs==0,
	S=find(~any(isnan(A),2));
    else
	S = [];
	% S = TissSpec(A(I,:),ncpgs);
	[S0,S,~,ind] = TissSpecBlocks(A(I,:),ncpgs,chrs,chr(I),pos(I),flank);
	% [S0,S,~,ind] = TissSpecBlocksResample(A(I,:),ncpgs,chrs,chr(I),pos(I),flank);
	% for i=1:10, SS=CompSens(A(I,:),S); S=[S;SS]; end;
	[S2,~,ind2] = ECC2(A(I,:),nECC,S,nECCmax); S = [S;S2];
	% S = randperm(length(I),1000);
	S = I(S);

	% count blocks
	s=0; for i=1:22, s=s+sum(diff(sort(pos(S(chr(S)==i))))>500); end;
	fprintf('%d CpGs ignored (due to NaNs) -- %d selected -- %d unique -- %d blocks\n', ...
		sum(any(isnan(A),2)), length(S), length(unique(S)), s)

	if 1,
	    % dump CpGs to file
	    fid=fopen(sprintf('CpGs.%dbp-block.%d.X1.xls',2*flank,ncpgs),'w');
	    fprintf(fid,'acc\tchr\tpos\tfrom\tto\tgroup\tname');
	    fprintf(fid,'\t%s',HD{:});
	    fprintf(fid,'\n');
	    for i=1:d,
		II=find(ind==i);
		for j=1:length(II),
		    s=S(II(j));
		    fprintf(fid,'%s\t%s\t%d\t%d\t%d\t%d\t%s',cg{s},chrs{chr(s)},pos(s),pos(s)-flank,pos(s)+flank, i, HD{i});
		    fprintf(fid,'\t%.1f%%', 100*A(s,:));
		    fprintf(fid,'\n');
		end
	    end

	    for i=1:d,
		II=find(ind==-i);
		for j=1:length(II),
		    s=S(II(j));
		fprintf(fid,'%s\t%s\t%d\t%d\t%d\t%d\t%s',cg{s},chrs{chr(s)},pos(s),pos(s)-flank,pos(s)+flank,-i,HD{i});
		    fprintf(fid,'\t%.1f%%', 100*A(s,:));
		    fprintf(fid,'\n');
		end
	    end
	    fclose(fid);

	    % dump CpGs to file
	    fid=fopen(sprintf('CpGs.ECC-pairs.%d.X1.xls',nECC),'w');
	    fprintf(fid,'acc\tchr\tpos\tfrom\tto\tname\t%%meth\tname\t%%meth'); fprintf(fid,'\n');
	    for i=1:length(S2),
		s=I(S2(i)); 
		fprintf(fid,'%s\t%s\t%d\t%d\t%d\t%s\t%.1f%%\t%s\t%.1f%%\n', ...
			cg{s}, chrs{chr(s)}, pos(s), pos(s)-flank, pos(s)+flank, ...
			HD{ind2(i,2)}, 100*A(s,ind2(i,2)), ...
			HD{ind2(i,3)}, 100*A(s,ind2(i,3)));
	    end
	    fclose(fid);
	end

    end
else
    disp('no feature selection');
end

if 1,
    figure(18); clf;
    imagesc(A(S(1:length(ind)),:),[0 1]);
    set(gca,'XTick',1:d,'XTickLabels',HD,'FontSize',12,'YTick',[1 5e3 10e3 15e3],'YTickLabels',{'0','5K','10K','15K'});
    XTickLabel_rotate;
    h=colorbar; set(h,'YLim',[0 1],'Ticks',[0 1],'TickLabels',{'0%','100%'},'FontSize',12);
    set(gcf,'PaperSize',[8 5],'PaperPosition',[0 0 8 5]); print(gcf,'-dpdf','-r300','Fig1A.pdf');

    figure(19); clf;
    imagesc(A(S(1:length(ind)),:),[0 1]); set(gca,'XTick',[],'YTick',[]);
    set(gcf,'PaperSize',[4 4],'PaperPosition',[0 0 4 4]); print(gcf,'-dpdf','-r300','Fig1B.pdf');

    % figure(20); clf;
    % imagesc(Data(S(1:length(ind)),J(1)),[0 1]); set(gca,'XTick',[],'YTick',[]);
    % set(gcf,'PaperSize',[.4 4],'PaperPosition',[0 0 .4 4]); print(gcf,'-dpdf','-r300','Fig1B2.pdf');
end
