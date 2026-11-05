# Existing aggregate dose contrasts, plotted without re-estimation.
suppressPackageStartupMessages(library(grid))
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
source <- normalizePath(args[1], winslash='/'); out <- args[2]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
d <- read.csv(source,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding='UTF-8-BOM')
metrics <- c('standardized_missed','false_alarm','cuff_conditioned')
rows <- c(53,48.5,41.8,37.3,23.5,19)
selected <- list();k<-0
for(m in metrics) for(db in c('mimic','eicu')) {
 k<-k+1
 ix <- which(d$metric==m & d$database==db & d$contrast=='highest_minus_lowest' & d$weighting=='pair' & d$standard_target==ifelse(m=='standardized_missed','pooled','not_applicable'))
 stopifnot(length(ix)==1)
 rr<-d[ix,];for(nm in c('estimate_pp','ci_lower_pp','ci_upper_pp'))rr[[nm]]<-as.numeric(rr[[nm]])
 stopifnot(rr$ci_lower_pp<=rr$estimate_pp,rr$ci_upper_pp>=rr$estimate_pp)
 selected[[k]]<-data.frame(element=paste(db,m,sep='_'),source_line=ix+1,estimate_pp=rr$estimate_pp,ci_lower_pp=rr$ci_lower_pp,ci_upper_pp=rr$ci_upper_pp,y_mm=rows[k])
}
sel<-do.call(rbind,selected)
fmt<-function(x) gsub('-', '\u2212',sprintf('%.1f',x),fixed=TRUE)
draw<-function(gray=FALSE) {
 grid.newpage();pushViewport(viewport(xscale=c(0,170),yscale=c(0,74)))
 tx<-function(label,x,y,size=7.5,col='#202020',just='left',bold=FALSE) grid.text(label,x,y,default.units='native',just=just,gp=gpar(fontfamily='Arial',fontsize=size,col=col,fontface=ifelse(bold,'bold','plain'),lineheight=1.25))
 ln<-function(x,y,col='#202020',lw=.8,lty=1) grid.lines(x,y,default.units='native',gp=gpar(col=col,lwd=lw/0.75,lty=lty,lineend='butt'))
 xp<-function(x) 60+(x+22)/54*68
 ln(c(60,128),c(14,14));ln(c(60,60),c(14,61))
 ln(rep(xp(0),2),c(14,61),col='#c2c2c2',lw=.7,lty='dashed')
 for(v in c(-20,-10,0,10,20,30)){ln(rep(xp(v),2),c(14,14-3.5*25.4/72));tx(gsub('-', '\u2212',as.character(v),fixed=TRUE),xp(v),10.2,just='centre')}
 tx('Highest minus lowest NE-eq stratum (percentage points)',94,4,just='centre')
 tx('Discrimination',3,58,size=8.1,bold=TRUE)
 tx('Meaning of an acceptable reading',3,29.8,size=8.1,bold=TRUE)
 for(y in c(56,27.5))ln(c(3,167),rep(y,2),col='#ededed',lw=.45)
 tx('Missed hypotension,\ndepth-standardized',4,50.75)
 tx('False alarm',4,39.55)
 tx('Arterial MAP <65 mmHg\ngiven cuff MAP \u226565 mmHg',4,21.25)
 tx('Difference (95% CI)',132,58)
 tx('\u2190 Lower at highest dose',xp(0)-1,64,size=6.4,col='#666666',just='right')
 tx('Higher at highest dose \u2192',xp(0)+1,64,size=6.4,col='#666666')
 for(j in 1:2){x<-c(128,153)[j];cl<-c('#202020','#787878')[j];grid.points(x,70,pch=c(16,15)[j],size=unit(1.35,'mm'),default.units='native',gp=gpar(col=cl,fill=cl));tx(c('MIMIC-IV','eICU')[j],x+2.5,70)}
 for(i in seq_len(nrow(sel))){
  r<-sel[i,];cl<-if(i%%2==1)'#202020' else '#787878'
  ln(xp(c(r$ci_lower_pp,r$ci_upper_pp)),rep(r$y_mm,2),col=cl,lw=.85)
  grid.points(xp(r$estimate_pp),r$y_mm,pch=if(i%%2==1)16 else 15,size=unit(1.35,'mm'),default.units='native',gp=gpar(col=cl,fill=cl))
  ss<-strsplit(fmt(r$estimate_pp),'.',fixed=TRUE)[[1]]
  tx(ss[1],136.5,r$y_mm,col=cl,just='right');tx(paste0('.',ss[2]),136.5,r$y_mm,col=cl)
  tx(paste0('(',fmt(r$ci_lower_pp),' to ',fmt(r$ci_upper_pp),')'),142,r$y_mm,col=cl)
 }
 popViewport()
}
svglite::svglite(file.path(out,'Figure2.svg'),width=170/25.4,height=74/25.4);draw();dev.off()
grDevices::cairo_pdf(file.path(out,'Figure2.pdf'),width=170/25.4,height=74/25.4,family='Arial');draw();dev.off()
ragg::agg_tiff(file.path(out,'Figure2.tif'),width=170,height=74,units='mm',res=600,compression='lzw');draw();dev.off()
ragg::agg_png(file.path(out,'Figure2.png'),width=170,height=74,units='mm',res=600);draw();dev.off()
ragg::agg_png(file.path(out,'Figure2_final_size.png'),width=170,height=74,units='mm',res=150);draw();dev.off()
ragg::agg_png(file.path(out,'Figure2_grayscale.png'),width=170,height=74,units='mm',res=150);draw(TRUE);dev.off()
tr<-do.call(rbind,lapply(seq_len(nrow(sel)),function(i)do.call(rbind,lapply(c('estimate_pp','ci_lower_pp','ci_upper_pp'),function(col)data.frame(figure='Figure2',element=sel$element[i],source_file=source,source_line=sel$source_line[i],source_column=col,source_value=sel[[col]][i],displayed_value=fmt(sel[[col]][i]))))))
write.csv(tr,file.path(out,'Figure2_number_trace.csv'),row.names=FALSE,fileEncoding='UTF-8')
cat('FIGURE2_EXPORT_PASS\n')
