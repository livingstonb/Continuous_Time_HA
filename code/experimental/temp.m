xlxpath = '/home/brian/Documents/temp/output_table.xlsx';
T = readtable(xlxpath);

a_lb = cellfun(@(x) str2num(convertStringsToChars(x)), table2cell(T(2,2:end)));
a_lb = a_lb(:);

chi0 = cellfun(@(x) str2num(convertStringsToChars(x)), table2cell(T(3,2:end)));
chi0 = chi0(:);

chi1 = cellfun(@(x) str2num(convertStringsToChars(x)), table2cell(T(4,2:end)));
chi1 = chi1(:);

chi2 = cellfun(@(x) str2num(convertStringsToChars(x)), table2cell(T(5,2:end)));
chi2 = chi2(:);

chivar = cellfun(@(x) str2num(convertStringsToChars(x)), table2cell(T(6,2:end)));
chivar = chivar(:);

ihtm = find(strcmp(table2cell(T(:,1)), 'L Wealth <= (1/6) Own Quart Inc'));
htm = cellfun(@(x) str2num(convertStringsToChars(x)), table2cell(T(ihtm,2:end)));
htm = htm(:);

iratio = find(strcmp(table2cell(T(:,1)), 'Wealthy HtM / Total HtM (at 1/6 qincome)'));
ratio = cellfun(@(x) str2num(convertStringsToChars(x)), table2cell(T(iratio,2:end)));
ratio = ratio(:);

imedian_liq = find(strcmp(table2cell(T(:,1)), 'Median Liq Assets'));
median_liq = cellfun(@(x) str2num(convertStringsToChars(x)), table2cell(T(imedian_liq,2:end)));
median_liq = median_liq(:);

imedian_tot = find(strcmp(table2cell(T(:,1)), 'Median Total Assets'));
median_tot = cellfun(@(x) str2num(convertStringsToChars(x)), table2cell(T(imedian_tot,2:end)));
median_tot = median_tot(:);

% Find indices that didn't converge
v = [median_liq - 0.1, median_tot - 1.7];
v = sum(v .^ 2, 2);
isuccess = v < 1e-8;

a_lb = a_lb(isuccess);
chi0 = chi0(isuccess);
chi1 = chi1(isuccess);
chi2 = chi2(isuccess);
chivar = chivar(isuccess);
htm = htm(isuccess);
ratio = ratio(isuccess);

table(a_lb, chi0, chi1, chi2, chivar, htm, ratio)