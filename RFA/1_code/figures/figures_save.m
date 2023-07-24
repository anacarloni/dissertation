hgsave([figname,'.fig'])
% print(figname,'-dpng')
% print(figname,'-depsc')
print(figname,'-dpdf','-bestfit')
% status = movefile(strcat(figname,'.eps'),strcat(dir,figname,'.eps'),'f');
% status = movefile(strcat(figname,'.png'),strcat(dir,figname,'.png'),'f');
status = movefile(strcat(figname,'.fig'),strcat(dir,figname,'.fig'),'f');
status = movefile(strcat(figname,'.pdf'),strcat(dir,figname,'.pdf'),'f');