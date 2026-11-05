# Redraw existing full-precision aggregates; no clinical analysis is performed.
# R-only drawing and export, using individually specified grid primitives.
# Usage: Rscript draw_figure1.R PROJECT_ROOT OUTPUT_ROOT PYTHON_FOR_GATE
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3)
invisible(Sys.setlocale("LC_CTYPE", "English_United States.UTF-8"))
root <- normalizePath(args[1],winslash="/",mustWork=TRUE)
out <- normalizePath(args[2],winslash="/",mustWork=FALSE)
python <- args[3]
figdir <- file.path(out,"figures"); review <- file.path(out,"review")
dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
dir.create(review,recursive=TRUE,showWarnings=FALSE)
library(grid)
for (pkg in c("svglite","ragg","jsonlite","systemfonts")) stopifnot(requireNamespace(pkg,quietly=TRUE))
stopifnot(any(systemfonts::system_fonts()$family=="Arial"))
sources <- c(s5="phase5/phase5_S5_counterfactual.csv",s3="phase3/revision4/data/S3_cuff_conditioned_probabilities.csv")
tables <- lapply(sources,function(s) read.csv(file.path(root,s),colClasses="character",check.names=FALSE,fileEncoding="UTF-8-BOM"))
strata <- c("<=0.1",">0.1-0.25",">0.25-0.5",">0.5-1.0",">1.0")
labels <- c("≤0.1",">0.1–0.25",">0.25–0.5",">0.5–1.0",">1.0")
bins <- c("65-69","70-74","75-79",">=80")
groups <- c("<=0.25",">0.25-1.0",">1.0")
trace <- list()
add <- function(panel,db,series,st,cat,key,match,val,lo=NULL,hi=NULL) {
  tab <- tables[[key]]; idx <- which(Reduce(`&`,lapply(names(match),function(k) tab[[k]]==match[[k]])))
  stopifnot(length(idx)==1); r <- tab[idx,,drop=FALSE]
  patients <- as.numeric(r[[if(key=="s5")"n_patients" else "denominator_patients"]])
  pairs <- as.numeric(r[[if(key=="s5")"n_pairs" else "denominator"]])
  stopifnot(patients>=10,pairs>=10)
  trace[[length(trace)+1]] <<- data.frame(panel=panel,database=db,series=series,stratum=st,category=cat,
    estimate=r[[val]],ci_lower=if(is.null(lo))"" else r[[lo]],ci_upper=if(is.null(hi))"" else r[[hi]],
    denominator_patients=patients,n_pairs=pairs,source_csv=sources[[key]],source_data_row=idx,source_line=idx+1,
    source_columns=paste(c(val,lo,hi),collapse=";"),stringsAsFactors=FALSE)
}
for (i in 1:2) {
  db <- c("mimic","eicu")[i]; name <- c("MIMIC-IV","eICU")[i]
  for (series in c("observed","expected","prevalence")) for (st in strata) {
    key <- if(series=="prevalence")"arterial_hypotension_prevalence" else series
    add(LETTERS[i],db,series,st,"","s5",list(database=db,analysis="primary",metric="art_lt65_given_cuff_ge65",stratum=st),
        key,if(series=="prevalence")NULL else paste0(key,"_ci_lower"),if(series=="prevalence")NULL else paste0(key,"_ci_upper"))
  }
  for(st in groups) for(cat in bins) add(LETTERS[i+2],db,st,st,cat,"s3",
    list(database=name,scenario="Primary five-minute pairing",stratum=paste0("dosegroup:",st),category=cat),
    "estimate","ci_lower","ci_upper")
}
data <- do.call(rbind,trace);stopifnot(nrow(data)==54)
datafile <- file.path(review,"Figure1_plot_data.csv")
write.csv(data,datafile,row.names=FALSE,fileEncoding="UTF-8",na="")
gate <- function(paths,stage) {
  log <- file.path(review,paste0("Figure1_gate_",stage,".log"))
  receipt <- file.path(review,paste0("Figure1_gate_",stage,".json"))
  status <- system2(python,c(shQuote(file.path(root,"tools/small_cell_gate.py")),shQuote(paths),"--report",shQuote(receipt)),stdout=log,stderr=log)
  stopifnot(status==0,jsonlite::fromJSON(receipt)$status=="PASS")
}
gate(datafile,"before")
contract <- list(conclusion="An acceptable single cuff reading carries greater residual arterial hypotension risk at higher vasopressor doses.",
  archetype="quantitative grid",backend="R grid",width_mm=170,height_mm=149.324324324324,
  font="Arial",body_pt=7.5,axes_linewidth_pt=.8,tick_length_pt=3.5,original_producer="phase3/icm_revision4/code/02_figure.py",
  source_csv=sources,analysis_reruns=0,source_values_unchanged=TRUE)
jsonlite::write_json(contract,file.path(review,"Figure1_contract.json"),pretty=TRUE,auto_unbox=TRUE)
width <- 170; height <- 149.324324324324
to_mm <- function(pt) pt*25.4/72
text <- function(s,x,y,size=7.5,colour="black",face="plain",just="centre",rot=0) grid.text(s,unit(x,"mm"),unit(y,"mm"),
  just=just,rot=rot,gp=gpar(fontfamily="Arial",fontsize=size,col=colour,fontface=face))
line <- function(x,y,colour="black",lw=.8,lt="solid") grid.lines(unit(x,"mm"),unit(y,"mm"),gp=gpar(col=colour,lwd=lw/.75,lty=lt,lineend="butt"))
point <- function(x,y,pch,col,fill=col,sz=3.6) grid.points(unit(x,"mm"),unit(y,"mm"),pch=pch,size=unit(sz,"pt"),gp=gpar(col=col,fill=fill,lwd=.8/.75))
cols <- c(observed="#0072B2",expected="#D55E00",prevalence="#8C8C8C","<=0.25"="#0072B2",">0.25-1.0"="#D55E00",">1.0"="#009E73")
pchs <- c(observed=21,expected=22,prevalence=23,"<=0.25"=21,">0.25-1.0"=22,">1.0"=24)
geometry <- list()
draw <- function(grayscale=FALSE,record=FALSE) {
  grid.newpage();pushViewport(viewport(width=unit(width,"mm"),height=unit(height,"mm")))
  usecols <- cols
  if(grayscale) for(k in names(usecols)) {v<-col2rgb(usecols[k])/255;usecols[k]<-gray(.2126*v[1]+.7152*v[2]+.0722*v[3])}
  for(row in 0:1) for(col in 0:1) {
    panel <- LETTERS[1+row*2+col]; left<-16+82*col; bottom<-if(row==0)94 else 25; pw<-65;ph<-43
    xmin<-if(row==0)-.4 else -.2;xmax<-if(row==0)4.4 else 3.2;ymax<-if(row==0)65 else 75
    xx<-function(x) left+(x-xmin)/(xmax-xmin)*pw; yy<-function(y) bottom+y/ymax*ph
    line(c(left,left,left+pw),c(bottom+ph,bottom,bottom),lw=.8)
    xt<-if(row==0)0:4 else 0:3; xl<-if(row==0)labels else c("65–69","70–74","75–79","≥80")
    for(j in seq_along(xt)) {line(rep(xx(xt[j]),2),c(bottom,bottom-to_mm(3.5)));text(xl[j],xx(xt[j]),bottom-4.2)}
    for(y in c(0,20,40,60)) {line(c(left,left-to_mm(3.5)),rep(yy(y),2));text(as.character(y),left-2.5,yy(y),just="right")}
    text(if(row==0)"NE-eq (µg/kg/min)" else "Cuff MAP (mmHg)",left+pw/2,bottom-9.2)
    if(col==0)text("Arterial MAP <65 mmHg (%)",4.3,bottom+ph/2,rot=90)
    if(row==0)text(c("MIMIC-IV","eICU")[col+1],left+pw/2,bottom+ph+6.5,size=8)
    text(panel,left,bottom+ph+6.5,size=8.5,face="bold",just="left")
    for(s in if(row==0)c("observed","expected","prevalence") else groups) {
      rr<-which(data$panel==panel & data$series==s);off<-if(s=="observed")-.08 else if(s=="expected").08 else 0
      xs<-seq_along(rr)-1+off;ys<-as.numeric(data$estimate[rr])*100
      cc<-usecols[s];line(xx(xs),yy(ys),cc,.85,if(s=="prevalence")"dotted" else "solid")
      for(j in seq_along(rr)) {
        n<-rr[j];px<-xx(xs[j]);py<-yy(ys[j])
        if(data$ci_lower[n]!="") {
          yl<-yy(as.numeric(data$ci_lower[n])*100);yh<-yy(as.numeric(data$ci_upper[n])*100)
          line(rep(px,2),c(yl,yh),cc,.8)
          line(c(px-to_mm(1.65),px+to_mm(1.65)),rep(yl,2),cc,.8)
          line(c(px-to_mm(1.65),px+to_mm(1.65)),rep(yh,2),cc,.8)
        } else {yl<-NA_real_;yh<-NA_real_}
        point(px,py,pchs[s],cc,if(s=="expected")"white" else cc,if(s=="prevalence")3.2 else 3.6)
        if(record)geometry[[length(geometry)+1]]<<-data.frame(trace_row=n,panel=panel,x_mm=px,y_mm=py,lower_y_mm=yl,upper_y_mm=yh,
           panel_left_mm=left,panel_bottom_mm=bottom,panel_width_mm=pw,panel_height_mm=ph,ymax=ymax)
      }
    }
  }
  lg <- function(s,label,x,y) {cc<-usecols[s];line(c(x,x+5),c(y,y),cc,.85,if(s=="prevalence")"dotted" else "solid");point(x+2.5,y,pchs[s],cc,if(s=="expected")"white" else cc);text(label,x+7,y,just="left")}
  lg("observed","Observed",20,77)
  lg("expected","Common-relation expected",55,77)
  lg("prevalence","Arterial hypotension prevalence",112,77)
  text("NE-eq (µg/kg/min)",49,7.2,just="right")
  lg("<=0.25","≤0.25",55,7.2)
  lg(">0.25-1.0",">0.25–1.0",85,7.2)
  lg(">1.0",">1.0",127,7.2)
  popViewport()
}
svglite::svglite(file.path(figdir,"Figure1.svg"),width=width/25.4,height=height/25.4,bg="white");draw(record=TRUE);dev.off()
grDevices::cairo_pdf(file.path(figdir,"Figure1.pdf"),width=width/25.4,height=height/25.4,family="Arial",bg="white");draw();dev.off()
ragg::agg_tiff(file.path(figdir,"Figure1.tif"),width=width,height=height,units="mm",res=600,compression="lzw",background="white");draw();dev.off()
ragg::agg_png(file.path(figdir,"Figure1.png"),width=width,height=height,units="mm",res=600,background="white");draw();dev.off()
for(res in c(150,96)) {
  file <- if(res==150)"Figure1_preview_170mm.png" else "Figure1_preview_100percent.png"
  ragg::agg_png(file.path(review,file),width=width,height=height,units="mm",res=res,background="white");draw();dev.off()
}
ragg::agg_png(file.path(review,"Figure1_grayscale_170mm.png"),width=width,height=height,units="mm",res=150,background="white");draw(grayscale=TRUE);dev.off()
write.csv(do.call(rbind,geometry),file.path(review,"Figure1_geometry.csv"),row.names=FALSE,fileEncoding="UTF-8",na="")
gate(c(file.path(figdir,paste0("Figure1.",c("pdf","svg","png","tif"))),datafile),"after")
cat(jsonlite::toJSON(list(status="PASS",records=nrow(data),width_mm=width,height_mm=height,font="Arial",gate_before="PASS",gate_after="PASS",R1_exceptions=0),auto_unbox=TRUE))
