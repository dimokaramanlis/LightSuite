function [cout, sigout, volout, stout, nameout] = ...
    reorganizeAreas(counts, signals, volumes, st, aggtype)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------
[Ngroups, ~, Nmice] = size(counts);
smat    = atlasTreeToMat(st);
smat    = smat(1:Ngroups, :);

if aggtype == 3
    agglist = st.id(1:Ngroups);
else
    agglist       = getAggList(aggtype); % based on Andy Peters selection
end
agglistcoarse = getAggList(1); % based on Andy Peters selection
coarsenames   = st.name(ismember(st.id, agglistcoarse));
%--------------------------------------------------------------------------
% we first organize areas in a reasonable depth

cout   = nan(numel(agglist), Nmice);
sigout  = nan(numel(agglist), Nmice);
volout  = nan(numel(agglist), 1);
nameout    = cell(numel(agglist), 3);
smatfin     = nan(numel(agglist), size(smat, 2));
for ii = 1:numel(agglist)

    [row, col]  = find(smat== agglist(ii));
    % counts are averaged in hemispheres
    countsignal = 2*squeeze(mean(counts(row, :, :), 2, 'omitmissing'));
    % countsignal = 2*squeeze(max(counts(row, :, :), [], 2, 'omitmissing'));

    countsignal = reshape(countsignal, numel(row), Nmice);
    cout(ii, :) = sum(countsignal, 1, 'omitmissing');
    % volumes are the whole volme
    volssel = volumes(row);
    volout(ii, :)   = sum(volssel, 1, 'omitmissing');
    
    % signal is averaged based on volume area
    wtscurr   = volssel/sum(volssel);
    bkgsignal = squeeze(mean(signals(row, :, :), 2, 'omitmissing'));
    bkgsignal = reshape(bkgsignal, numel(row), Nmice);
    signaluse = sum(wtscurr.*bkgsignal, 1, 'omitmissing');
    sigout(ii,:) = signaluse;

    idsgroup   = find(ismember(st.id, st.id(row)));
    structname = st.name{idsgroup(1)};
    structabbr = st.acronym{idsgroup(1)};
    
    idcoarse = ismember(agglistcoarse, smat(idsgroup(1),:));
    if nnz(idcoarse)>0
        lastname = coarsenames{idcoarse};
    else
        lastname = '';
    end
    % structname = st.name(st.id == mode(st.parent_structure_id(row)));
    % structabbr = st.acronym(st.id == mode(st.parent_structure_id(row)));
    nameout(ii, :) = {structname, structabbr, lastname};
    smatfin(ii,: )  = smat(idsgroup(1),:);

end

stout = st(ismember(st.id, agglist),:);
%--------------------------------------------------------------------------
% we then get rid of things 
% 
% 1. without volume in the atlas
% 2. Hindbrain, Cerebellum, Olfactory bulb, fiber tracts

iclassrem = contains(lower(st.name), 'cerebellum') |...
    contains(lower(st.name), 'fiber tracts') |...
    contains(lower(st.name), 'hindbrain') | ...
    contains(lower(st.name), 'olfactory bulb')|...
    contains(lower(st.name), 'ventricle');

ikeeptype   = ~any(ismember(smat, st.id(find(iclassrem))), 2);
ikeepvolume = volout > 0;

ikeep      = ikeepvolume; % & ikeeptype;
cout  = cout(ikeep,  :);
volout = volout(ikeep, :);
sigout = sigout(ikeep, :);
nameout = nameout(ikeep, :);
stout      = stout(ikeep, :);
% smatnew    = smat(ikeep, :);

%--------------------------------------------------------------------------

% areaids   = sorg(toindex);
% [unids, ~, ic] = unique(areaids);
% Nareas         = numel(unids);
% 
% finalvolumes = accumarray(ic, volumesnew, [Nareas 1], @sum);
% finalcounts  = nan(Nareas, 2, Nmice, 'single');
% for imouse = 1:Nmice
%     for ii = 1:2 
%         finalcounts(:, ii, imouse) = accumarray(ic, countsnew(:, ii, imouse), [Nareas 1], @nansum);
%     end
% end
% finalids = lower(st.name(findfirst(unids == st.id', 2)));
% %--------------------------------------------------------------------------
% 
% ivis = contains(finalids,'primary motor')|...
%     contains(finalids,'primary visual')|...
%     contains(finalids,'prelimbic')|...
%     contains(finalids,'infralimbic')|...
%     contains(finalids,'primary auditory')|...
%     contains(finalids,'secondary motor')|...
%     contains(finalids,'cingulate');
% 
% 
% tomerge  = find(~isnan(smat(:, depth + 1)));
% lowdepth = find(isnan(smat(:, depth + 1)));
% %--------------------------------------------------------------------------
% % we first keep the lowdepth areas
% countslowdepth  = counts(lowdepth(volumes(lowdepth) > 0), :, :);
% volumeslowdepth = volumes(lowdepth(volumes(lowdepth) > 0));
% 
% 
% %--------------------------------------------------------------------------
% 
% smat()
% 
% [newcheck,~,ic] = unique(smat(:,depth));
% ikeep   = ~isnan(newcheck);
% newvols = accumarray(ic, volumes, [], @sum);
% newvols = newvols(ikeep);
% 
% locations = findIndicesLocation(newcheck(ikeep), st.id(1:Ngroups));
% 
% st.name(locations)
% 
% ismember(st.id(1:Ngroups), newcheck(ikeep))
% 

end
    
function aggList = getAggList(aggType)

switch aggType
    case 1 % 10 top-level divisions: {'Isocortex', 'HPF', 'OLF', 'CTXsp', 'CNU', 'TH', 'HY', 'MB', 'HB', 'CB'}
        % aggList = [315 1089 698 703 623 549 1097 313 1065 512];
        %10 top-level divisions: {'Isocortex', 'HPF', 'OLF', 'CTXsp', 'CNU', 'TH', 'HY', 'MB', 'HB', 'CB'}
        aggList = [315  698 1089 703 623 549 1097 313];
        
    case 2 % 316 manually selected regions
        % aggList = [184 985 993 353 329 337 345 369 361 182305689 ...
        %             378 1057 677 1011 1002 1027 1018 402 394 ...
        %             409 385 425 533 312782574 312782628 39 48 972 44 723 ...
        %             731 738 746 104 111 119 894 879 886 312782546 417 ...
        %             541 922 895 507 151 159 597 605 814 961 619 639 647 ...
        %             788 566 382 423 463 726 982 19 918 926 843 1037 1084 ...
        %             502 484682470 589508447 484682508 583 952 966 131 295 ...
        %             319 780 672 56 998 754 250 258 266 310 333 23 292 536 ...
        %             1105 403 1022 1031 342 298 564 596 581 351 629 685 718 ...
        %             725 733 741 563807435 406 609 1044 475 170 218 1020 1029 ...
        %             325 560581551 255 127 64 1120 1113 155 59 362 366 1077 ...
        %             149 15 181 560581559 189 599 907 575 930 560581563 262 ...
        %             27 563807439 178 321 483 186 390 38 30 118 ...
        %             223 72 263 272 830 452 523 1109 126 133 347 286 338 ...
        %             576073699 689 88 210 491 525 557 515 980 1004 63 693 ...
        %             946 194 226 364 576073704 173 470 614 797 302 4 580 271 ...
        %             874 381 749 607344830 246 128 294 795 215 531 ...
        %             628 634 706 1061 549009203 616 214 35 549009211 975 115 ...
        %             606826663 757 231 66 75 58 374 1052 12 100 197 591 872 ...
        %             612 7 867 398 280 880 599626927 898 931 1093 318 534 574 ...
        %             621 549009215 549009219 549009223 549009227 679 147 162 ...
        %             604 146 238 350 358 207 96 101 711 1039 903 642 651 429 ...
        %             437 445 589508451 653 661 135 839 1048 372 83 136 106 ...
        %             203 235 307 395 852 859 938 177 169 995 1069 209 202 225 ...
        %             217 765 773 781 206 230 222 912 976 984 1091 936 944 951 ...
        %             957 968 1007 1056 1064 1025 1033 1041 1049 989 91 846 ...
        %             589508455];

         aggList = [184 985 993 353 329 337 345 369 361 182305689 ...
                    378 1057 677 1011 1002 1027 1018 402 394 ...
                    409 385 425 533 312782574 312782628 39 48 972 44 723 ...
                    731 738 746 104 111 119 894 879 886 312782546 417 ...
                    541 922 895 507 151 159 597 605 814 961 619 639 647 ...
                    788 566 382 423 463 726 982 19 918 926 843 1037 1084 ...
                    502 484682470 589508447 484682508 583 952 966 131 295 ...
                    319 780 672 56 998 754 250 258 266 310 333 23 292 536 ...
                    1105 403 1022 1031 342 298 564 596 581 351 629 685 718 ...
                    725 733 741 563807435 406 609 1044 475 170 218 1020 1029 ...
                    325 560581551 255 127 64 1120 1113 155 59 362 366 1077 ...
                    149 15 181 560581559 189 599 907 575 930 560581563 262 ...
                    27 563807439 178 321 483 186 390 38 30 118 ...
                    223 72 263 272 830 452 523 1109 126 133 347 286 338 ...
                    576073699 689 88 210 491 525 557 515 980 1004 63 693 ...
                    946 194 226 364 576073704 173 470 614 797 302 4 580 271 ...
                    874 381 749 607344830 246 128 294 795 215 531 ...
                    628 634 706 1061 549009203 616 214 35 549009211 975 115 ...
                    606826663 757 231 66 75 58 374 1052 12 100 197 591 872 ];

end

end