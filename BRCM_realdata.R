.envnum <- function(k,d){v<-Sys.getenv(k);if(nzchar(v))as.numeric(v) else d}
.envchr <- function(k,d){v<-Sys.getenv(k);if(nzchar(v))v else d}
N_REPEATS <- .envnum("BRCM_REPEATS", 10)
N_FOLDS   <- .envnum("BRCM_FOLDS", 5)
OOF_K     <- .envnum("BRCM_OOFK", 3)
BASE_SEED <- .envnum("BRCM_SEED", 2000)
N_CORES   <- .envnum("BRCM_CORES", max(1, parallel::detectCores() - 1))
OUT_DIR   <- .envchr("BRCM_OUT", "BRCMreal")
PROG_CHUNKS <- as.integer(.envnum("BRCM_PROGRESS", 20))   
R2_THRESH <- .envnum("BRCM_R2THRESH", 0.01)
AUTORUN   <- toupper(.envchr("BRCM_AUTORUN", "TRUE")) == "TRUE"
STRICT_PKGS <- toupper(.envchr("BRCM_STRICT", "TRUE")) == "TRUE"
REQUIRED_PKGS  <- c("mlbench","mgcv","rpart","randomForest","e1071")
MANUSCRIPT_PKGS <- c("xgboost","lightgbm","catboost")          
has <- function(p) requireNamespace(p, quietly = TRUE)
.check_pkgs <- function(){
  need <- if (STRICT_PKGS) c(REQUIRED_PKGS, MANUSCRIPT_PKGS) else REQUIRED_PKGS
  miss <- need[!vapply(need, requireNamespace, TRUE, quietly = TRUE)]
  if (length(miss)) {
    msg <- paste0("Missing packages required for this run: ", paste(miss, collapse = ", "))
    if (STRICT_PKGS) stop(msg, ". Install them for the manuscript run, or set BRCM_STRICT=FALSE for a development run that skips missing methods.")
    warning(msg, " -- those methods will be skipped in this development run.")}
  invisible(TRUE)}

bern <- function(u,m){k<-0:m;u<-pmin(pmax(u,1e-9),1-1e-9);sapply(k,function(kk)dbinom(kk,m,u))}
u_train <- function(x){(rank(x,ties.method="average")-0.5)/length(x)}
u_test  <- function(xn,xt){                                
  xs<-sort(xt);n<-length(xs)
  n_le<-findInterval(xn,xs)                                
  n_lt<-findInterval(xn,xs,left.open=TRUE)                
  (n_lt+n_le)/(2*n)}                                      
.strat_folds <- function(y,K=3,seed=2026){
    old<-if(exists(".Random.seed",envir=.GlobalEnv))get(".Random.seed",envir=.GlobalEnv) else NULL
  set.seed(seed);f<-integer(length(y))
  for(cl in unique(y)){i<-which(y==cl);i<-i[sample.int(length(i))];f[i]<-rep(seq_len(K),length.out=length(i))}
  if(!is.null(old))assign(".Random.seed",old,envir=.GlobalEnv) else if(exists(".Random.seed",envir=.GlobalEnv))rm(".Random.seed",envir=.GlobalEnv)
  f}
.irls <- function(X,y,pen,maxit=100,tol=1e-8){
  p<-ncol(X);b<-rep(0,p);Pen<-diag(pen,p,p);conv<-FALSE;it<-0L
  for(it in 1:maxit){eta<-drop(X%*%b);mu<-1/(1+exp(-eta));w<-pmax(mu*(1-mu),1e-6)
    z<-eta+(y-mu)/w;A<-crossprod(X*sqrt(w))+Pen;g<-crossprod(X,w*z)
    bn<-tryCatch(solve(A,g),error=function(e)tryCatch(solve(A+1e-6*diag(p),g),error=function(e2)NULL))
    if(is.null(bn))break
    if(max(abs(bn-b))<tol){b<-bn;conv<-TRUE;break};b<-bn}
  list(beta=as.vector(b),converged=conv,iter=it)}
.ll <- function(y,p){p<-pmin(pmax(p,1e-12),1-1e-12);-mean(y*log(p)+(1-y)*log(1-p))}
.spw <- function(X,y,c_spear){vapply(seq_len(ncol(X)),function(j){
  r<-abs(suppressWarnings(stats::cor(X[,j],y,method="spearman")));if(is.na(r))r<-0;1/(r+c_spear)},numeric(1))}

build_main <- function(X,cat_mask,m,mode="bern"){
  cols<-list();meta<-list();blk<-c()
  for(j in seq_len(ncol(X))){
    if(cat_mask[j]){
      lv<-sort(unique(X[,j]))
      if(length(lv)<2L){meta[[j]]<-list(type="cat",levels=numeric(0),ctr=numeric(0),drop=TRUE);next}
      oth<-lv[-1]
      M<-vapply(oth,function(c)as.numeric(X[,j]==c),numeric(nrow(X)));if(is.null(dim(M)))M<-matrix(M,ncol=1)
      ctr<-colMeans(M);M<-sweep(M,2,ctr)
      meta[[j]]<-list(type="cat",ref=lv[1],levels=oth,ctr=ctr,drop=FALSE)
      cols[[length(cols)+1]]<-M;blk<-c(blk,rep(j,ncol(M)))
    } else if(mode=="raw"){
      mu<-mean(X[,j]);sdv<-stats::sd(X[,j])+1e-12
      meta[[j]]<-list(type="raw",mu=mu,sd=sdv,drop=FALSE)
      cols[[length(cols)+1]]<-matrix((X[,j]-mu)/sdv,ncol=1);blk<-c(blk,j)
    } else if(mode=="ranklin"){
      xt<-X[,j];u<-u_train(xt);um<-mean(u)
      meta[[j]]<-list(type="ranklin",xsort=sort(xt),umean=um,drop=FALSE)
      cols[[length(cols)+1]]<-matrix(u-um,ncol=1);blk<-c(blk,j)
    } else {
      xt<-X[,j];B<-bern(u_train(xt),m);ctr<-colMeans(B);B<-sweep(B,2,ctr)
      meta[[j]]<-list(type="cont",xsort=sort(xt),ctr=ctr,drop=FALSE)
      cols[[length(cols)+1]]<-B;blk<-c(blk,rep(j,ncol(B)))}}
  list(D=if(length(cols))do.call(cbind,cols) else matrix(0,nrow(X),0),block=blk,meta=meta)}

map_main_test <- function(Xte,meta,m){
  cols<-list()
  for(j in seq_along(meta)){mj<-meta[[j]];if(isTRUE(mj$drop))next
    if(mj$type=="cat"){
      M<-vapply(mj$levels,function(c)as.numeric(Xte[,j]==c),numeric(nrow(Xte)));if(is.null(dim(M)))M<-matrix(M,ncol=1)
      cols[[length(cols)+1]]<-sweep(M,2,mj$ctr)
    } else if(mj$type=="raw"){
      cols[[length(cols)+1]]<-matrix((Xte[,j]-mj$mu)/mj$sd,ncol=1)
    } else if(mj$type=="ranklin"){
      cols[[length(cols)+1]]<-matrix(u_test(Xte[,j],mj$xsort)-mj$umean,ncol=1)
    } else cols[[length(cols)+1]]<-sweep(bern(u_test(Xte[,j],mj$xsort),m),2,mj$ctr)}
  if(length(cols))do.call(cbind,cols) else matrix(0,nrow(Xte),0)}

.tensor <- function(uj,uk,m2){Bj<-bern(uj,m2);Bk<-bern(uk,m2)
  do.call(cbind,lapply(seq_len(ncol(Bj)),function(a)Bj[,a]*Bk))}
.pair_A_new <- function(Dnew,block,j,k) cbind(1,Dnew[,block==j,drop=FALSE],Dnew[,block==k,drop=FALSE])
.pair_T_new <- function(Xnew,meta,j,k,m2) .tensor(u_test(Xnew[,j],meta[[j]]$xsort),u_test(Xnew[,k],meta[[k]]$xsort),m2)
pair_pure <- function(j,k,X,meta,m,m2,D,block){
  A<-.pair_A_new(D,block,j,k)
  Tt<-.tensor(u_test(X[,j],meta[[j]]$xsort),u_test(X[,k],meta[[k]]$xsort),m2)
  Theta<-qr.coef(qr(A),Tt);Theta[is.na(Theta)]<-0
  list(T=Tt-A%*%Theta,Theta=Theta,j=j,k=k,ncol=ncol(Tt))}

brcm_fit <- function(X,y,cat_mask,m=10,m2=3,Kmax=4,lam_int=1,c_spear=0.05,
                     lam_grid=c(0.25,0.5,1,2,4),r2_thresh=R2_THRESH,min_n_inter=300,cv_K=3,top_feat=8,
                     screen_mode=c("nested","outer"),
                     use_rank=TRUE,use_bernstein=TRUE,use_spearman=TRUE,use_inter=TRUE){
  screen_mode<-match.arg(screen_mode)
  X<-as.matrix(X);n<-nrow(X);folds<-.strat_folds(y,cv_K)
  mode<-if(!use_rank)"raw" else if(!use_bernstein)"ranklin" else "bern"
  mm<-build_main(X,cat_mask,m,mode)
  n_valid<-sum(vapply(1:cv_K,function(f)length(unique(y[which(folds!=f)]))>=2,TRUE))
  cvll<-vapply(lam_grid,function(lam){e<-c();ok<-TRUE
    for(f in 1:cv_K){tri<-which(folds!=f);vai<-which(folds==f)
      if(length(unique(y[tri]))<2)next                        
      wf<-if(use_spearman) .spw(X[tri,,drop=FALSE],y[tri],c_spear) else rep(1,ncol(X))
      mf<-build_main(X[tri,,drop=FALSE],cat_mask,m,mode)
      if(ncol(mf$D)==0){ok<-FALSE;break}
      ft<-.irls(cbind(1,mf$D),y[tri],c(0,lam*wf[mf$block]));if(!ft$converged){ok<-FALSE;break}
      Dv<-map_main_test(X[vai,,drop=FALSE],mf$meta,m)
      e<-c(e,.ll(y[vai],1/(1+exp(-drop(cbind(1,Dv)%*%ft$beta)))))}
    if(ok&&length(e)==n_valid&&n_valid>0)mean(e) else NA_real_},numeric(1))
  lam<-if(all(is.na(cvll)))1 else lam_grid[which.min(cvll)]
  w<-if(use_spearman) .spw(X,y,c_spear) else rep(1,ncol(X))
  D<-mm$D;penblk<-w[mm$block]
  sel<-list();selinfo<-data.frame();scores<-numeric(0)
  cont_all<-which(!cat_mask & !vapply(mm$meta,function(z)isTRUE(z$drop),TRUE))
  .screen <- function(idx,rows){                      
    if(length(idx)<=top_feat) return(idx)
    rho<-vapply(idx,function(j){r<-abs(suppressWarnings(stats::cor(X[rows,j],y[rows],method="spearman")))
      if(is.na(r))0 else r},numeric(1))
    idx[order(rho,decreasing=TRUE)[seq_len(top_feat)]]}
  .pcan <- function(pr) as.integer(sort(as.integer(pr)))
  .pkey <- function(pr){q<-.pcan(pr);paste0(q[1],"_",q[2])}
  if(screen_mode=="outer"){
    cont<-.screen(cont_all,seq_len(n))
    pairs<-if(use_inter && use_rank && length(cont)>=2)
      lapply(combn(cont,2,simplify=FALSE),.pcan) else list()
    fold_ok<-NULL                                    
  } else {
    fold_pairs<-vector("list",cv_K)
    for(f in 1:cv_K){tri<-which(folds!=f)
      cf<-if(length(unique(y[tri]))<2) integer(0) else .screen(cont_all,tri)
      fold_pairs[[f]]<-if(use_inter && use_rank && length(cf)>=2)
        lapply(combn(cf,2,simplify=FALSE),.pcan) else list()}
    uni<-unique(unlist(lapply(fold_pairs,function(L)vapply(L,.pkey,"")),use.names=FALSE))
    pairs<-lapply(uni,function(k)as.integer(strsplit(k,"_",fixed=TRUE)[[1]]))
    fold_ok<-lapply(pairs,function(pr)vapply(fold_pairs,function(L)
      .pkey(pr) %in% vapply(L,.pkey,""),TRUE))
    names(fold_ok)<-vapply(pairs,.pkey,"")
  }
  n_screened_pairs<-length(pairs)
  interaction_search_enabled<-(length(pairs)>0L && n>=min_n_inter)
  if(interaction_search_enabled){
    scores<-vapply(seq_along(pairs),function(pi_){pr<-pairs[[pi_]]
      rv<-c();rh<-c();nf<-0L
      for(f in 1:cv_K){tri<-which(folds!=f);vai<-which(folds==f)
        if(length(unique(y[tri]))<2)next
        if(!is.null(fold_ok) && !isTRUE(fold_ok[[pi_]][f]))next
        nf<-nf+1L
        wf<-if(use_spearman) .spw(X[tri,,drop=FALSE],y[tri],c_spear) else rep(1,ncol(X))
        mf<-build_main(X[tri,,drop=FALSE],cat_mask,m);if(ncol(mf$D)==0)next
        ft<-.irls(cbind(1,mf$D),y[tri],c(0,lam*wf[mf$block]));if(!ft$converged)next
        r_tr<-y[tri]-1/(1+exp(-drop(cbind(1,mf$D)%*%ft$beta)))
        pp<-tryCatch(pair_pure(pr[1],pr[2],X[tri,,drop=FALSE],mf$meta,m,m2,mf$D,mf$block),error=function(e)NULL)
        if(is.null(pp))next
        g<-tryCatch(solve(crossprod(pp$T)+lam_int*diag(ncol(pp$T)),crossprod(pp$T,r_tr)),error=function(e)NULL)
        if(is.null(g))next
        Dv<-map_main_test(X[vai,,drop=FALSE],mf$meta,m)
        r_va<-y[vai]-1/(1+exp(-drop(cbind(1,Dv)%*%ft$beta)))
        Tv<-.pair_T_new(X[vai,,drop=FALSE],mf$meta,pr[1],pr[2],m2)
        Tvp<-Tv-.pair_A_new(Dv,mf$block,pr[1],pr[2])%*%pp$Theta
        rv<-c(rv,r_va);rh<-c(rh,drop(Tvp%*%g))}
      if(!is.null(fold_ok) && nf < (floor(cv_K/2)+1L)) return(-Inf)
      if(length(rv)<10)return(-Inf);1-sum((rv-rh)^2)/sum((rv-mean(rv))^2)},numeric(1))
    keep<-which(scores>r2_thresh);keep<-keep[order(scores[keep],decreasing=TRUE)];keep<-head(keep,Kmax)
    sel<-pairs[keep]
    if(length(keep))selinfo<-data.frame(j=vapply(sel,`[`,0L,1),k=vapply(sel,`[`,0L,2),
      oof_R2=round(scores[keep],5),
      n_folds_screened=if(is.null(fold_ok)) cv_K else vapply(keep,function(i)sum(fold_ok[[i]]),0L))}
  Dfull<-D;ipen<-c();inter<-list()
  for(pr in sel){pp<-tryCatch(pair_pure(pr[1],pr[2],X,mm$meta,m,m2,D,mm$block),error=function(e)NULL)
    if(is.null(pp))next
    Dfull<-cbind(Dfull,pp$T);ipen<-c(ipen,rep(lam_int,pp$ncol));inter[[length(inter)+1]]<-pp}
  ff<-.irls(cbind(1,Dfull),y,c(0,lam*penblk,ipen))
  structure(list(beta=ff$beta,lambda=lam,meta=mm$meta,block=mm$block,m=m,m2=m2,use_rank=use_rank,
    inter=inter,selinfo=selinfo,all_oof_R2=scores,screen_mode=screen_mode,
    n_screened_pairs=n_screened_pairs,interaction_search_enabled=interaction_search_enabled,
    cat_mask=cat_mask,mode=mode,converged=ff$converged,
    iter=ff$iter,cv_logloss=if(all(is.na(cvll)))NA_real_ else min(cvll,na.rm=TRUE),
    n_inter=length(inter)),class="brcm")}


brcm_components <- function(fit,Xnew){
  Xnew<-as.matrix(Xnew)
  D<-map_main_test(Xnew,fit$meta,fit$m)
  b<-fit$beta;M<-length(fit$block);main<-matrix(0,nrow(Xnew),length(fit$meta))
  for(j in unique(fit$block)){i<-which(fit$block==j);main[,j]<-D[,i,drop=FALSE]%*%b[1+i]}
  ints<-NULL;off<-1L+M
  if(length(fit$inter)){ints<-matrix(0,nrow(Xnew),length(fit$inter))
    for(q in seq_along(fit$inter)){it<-fit$inter[[q]]
      Tp<-.pair_T_new(Xnew,fit$meta,it$j,it$k,fit$m2)-.pair_A_new(D,fit$block,it$j,it$k)%*%it$Theta
      ints[,q]<-Tp%*%b[off+seq_len(it$ncol)];off<-off+it$ncol}}
  list(intercept=b[1],main=main,inter=ints)}
brcm_predict <- function(fit,Xnew){cp<-brcm_components(fit,Xnew)
  1/(1+exp(-(cp$intercept+rowSums(cp$main)+if(is.null(cp$inter))0 else rowSums(cp$inter))))}

minus_labels <- function(x) sub("^-", "\u2212", format(x, trim = TRUE))
axis_minus <- function(side, at = NULL, ...) {
  if (is.null(at)) at <- axTicks(side)
  axis(side, at = at, labels = minus_labels(at), ...)
}

brcm_plot_shapes <- function(fit,feature_names=NULL,file=NULL,ncol_grid=2){
  meta<-fit$meta;p<-length(meta);b<-fit$beta;blk<-fit$block
  if(is.null(feature_names))feature_names<-paste0("V",seq_len(p))
  if(!is.null(file))png(file,width=3000,height=2200,res=300)
  op<-par(mfrow=c(ceiling(p/ncol_grid),ncol_grid),mar=c(4.2,4.4,2.4,1));on.exit(par(op),add=TRUE)
  for(j in seq_len(p)){mj<-meta[[j]];i<-which(blk==j)
    if(isTRUE(mj$drop)||!length(i)){plot.new();title(main=paste0(feature_names[j]," (dropped)"));next}
    bj<-b[1+i]
    if(mj$type=="cat"){
      contrib<-c(0,bj)-sum(mj$ctr*bj)                     
      barplot(contrib,names.arg=c(as.character(mj$ref),as.character(mj$levels)),
        main=feature_names[j],ylab="Logit contribution",xlab="Category",col="#4C72B0",border=NA,yaxt="n")
      axis_minus(2)
      abline(h=0,lty=3,col="grey50")
    } else if(mj$type=="cont"){
      ug<-seq(0.005,0.995,length.out=200)
      plot(ug,drop(sweep(bern(ug,fit$m),2,mj$ctr)%*%bj),type="l",lwd=3,col="#C44E52",
        xlab="Normalized rank (ECDF)",ylab="Logit contribution",main=feature_names[j],yaxt="n")
      axis_minus(2)
      abline(h=0,lty=3,col="grey60")
    } else {plot.new();title(main=paste0(feature_names[j]," (linear)"))}}
  if(!is.null(file))dev.off();invisible(NULL)}

brcm_plot_interaction <- function(fit,which=1,feature_names=NULL,file=NULL,ngrid=40){
  if(length(fit$inter)<which){message("no such interaction");return(invisible(NULL))}
  it<-fit$inter[[which]];j<-it$j;k<-it$k
  if(is.null(feature_names))feature_names<-paste0("V",seq_len(length(fit$meta)))
  off<-1L+length(fit$block);if(which>1)for(q in 1:(which-1))off<-off+fit$inter[[q]]$ncol
  g<-fit$beta[off+seq_len(it$ncol)]                        
  ug<-seq(0.005,0.995,length.out=ngrid);gr<-expand.grid(uj=ug,uk=ug)
  Tp<-.tensor(gr$uj,gr$uk,fit$m2)-
      cbind(1,sweep(bern(gr$uj,fit$m),2,fit$meta[[j]]$ctr),
              sweep(bern(gr$uk,fit$m),2,fit$meta[[k]]$ctr))%*%it$Theta
  z<-matrix(drop(Tp%*%g),ngrid,ngrid)
  if(!is.null(file))png(file,width=2100,height=1800,res=300)
  op<-par(mar=c(4.5,4.5,3,2));on.exit(par(op),add=TRUE)
  image(ug,ug,z,xlab=paste0(feature_names[j]," (rank)"),ylab=paste0(feature_names[k]," (rank)"),
    main=sprintf("Pure interaction surface: %s x %s",feature_names[j],feature_names[k]),
    col=hcl.colors(50,"Blue-Red"));contour(ug,ug,z,add=TRUE,col="grey20")
  if(!is.null(file))dev.off();invisible(NULL)}

.auc <- function(y,p){n1<-sum(y==1);n0<-sum(y==0);if(n1==0||n0==0)return(NA_real_)
  r<-rank(p);(sum(r[y==1])-n1*(n1+1)/2)/(n1*n0)}
.safe_mean   <- function(x){x<-x[is.finite(x)];if(!length(x))NA_real_ else mean(x)}
.safe_median <- function(x){x<-x[is.finite(x)];if(!length(x))NA_real_ else stats::median(x)}
.name_seed <- function(nm){z<-utf8ToInt(nm);sum(z*seq_along(z))%%100000L}

.pr_auc <- function(y,p){                                  
  if(sum(y)==0)return(NA_real_)
  P<-sum(y)
  ord<-order(p,decreasing=TRUE);ps<-p[ord];ys<-y[ord]
  g<-cumsum(c(TRUE,diff(ps)!=0))                           
  tp_grp<-tapply(ys,g,sum);n_grp<-tapply(rep(1,length(ys)),g,sum)
  cum_tp<-cumsum(tp_grp);cum_n<-cumsum(n_grp)
  prec<-as.numeric(cum_tp/cum_n);rec<-as.numeric(cum_tp/P)
  sum(prec*diff(c(0,rec)))}
.ece <- function(y,p,B=10){b<-cut(p,breaks=seq(0,1,length.out=B+1),include.lowest=TRUE)
  s<-0;for(lv in levels(b)){i<-which(b==lv);if(!length(i))next
    s<-s+length(i)/length(y)*abs(mean(y[i])-mean(p[i]))};s}
.ece_adaptive <- function(y,p,B=10){
  n<-length(p); if(n<2L) return(NA_real_)
  B<-max(2L,min(B,floor(n/5))); if(B<2L) return(NA_real_)
  q<-stats::quantile(p,probs=seq(0,1,length.out=B+1),names=FALSE,type=7)
  q[1]<--Inf; q[length(q)]<-Inf
  b<-cut(p,breaks=unique(q),include.lowest=TRUE)
  s<-0;for(lv in levels(b)){i<-which(b==lv);if(!length(i))next
    s<-s+length(i)/n*abs(mean(y[i])-mean(p[i]))};s}
.ici <- function(y,p){
  n<-length(p); if(n<10L||length(unique(y))<2L) return(NA_real_)
  pc<-pmin(pmax(p,1e-6),1-1e-6)
  if(stats::sd(pc) < 1e-8) return(abs(mean(y) - mean(pc)))
  fit<-suppressWarnings(tryCatch(stats::loess(y~pc,span=0.75,degree=1,family="gaussian",
                             control=stats::loess.control(surface="direct")),
                error=function(e)NULL))
  ph<-if(!is.null(fit)) suppressWarnings(tryCatch(as.vector(stats::predict(fit,newdata=data.frame(pc=pc))),
                                                  error=function(e)NULL)) else NULL
  if(is.null(ph)||length(ph)!=n||any(!is.finite(ph))||stats::sd(ph)<1e-12){
    lo<-log(pc/(1-pc))
    g<-suppressWarnings(tryCatch(stats::glm(y~lo,family=stats::binomial()),error=function(e)NULL))
    if(is.null(g)) return(abs(mean(y)-mean(pc)))
    ph<-suppressWarnings(as.vector(stats::predict(g,type="response")))}
  if(any(!is.finite(ph))) return(abs(mean(y)-mean(pc)))
  ph<-pmin(pmax(ph,0),1)
  mean(abs(ph-pc))}
.cal <- function(y,p){p<-pmin(pmax(p,1e-6),1-1e-6);lo<-log(p/(1-p))
  f<-tryCatch(stats::glm(y~lo,family=stats::binomial()),error=function(e)NULL)
  if(is.null(f)||any(is.na(stats::coef(f))))return(c(NA_real_,NA_real_))
  c(stats::coef(f)[2],stats::coef(f)[1])}
.thr_metrics <- function(y,p,thr){yh<-as.integer(p>=thr)
  tp<-sum(yh==1&y==1);tn<-sum(yh==0&y==0);fp<-sum(yh==1&y==0);fn<-sum(yh==0&y==1)
  se<-if(tp+fn>0)tp/(tp+fn) else NA;sp<-if(tn+fp>0)tn/(tn+fp) else NA
  pr<-if(tp+fp>0)tp/(tp+fp) else NA
  f1<-if(!is.na(pr)&&!is.na(se)&&(pr+se)>0)2*pr*se/(pr+se) else NA
  den<-sqrt(as.numeric(tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  mcc<-if(den>0)(as.numeric(tp)*tn-as.numeric(fp)*fn)/den else NA
  acc<-(tp+tn)/length(y);po<-acc
  pe<-((tp+fp)*(tp+fn)+(tn+fn)*(tn+fp))/length(y)^2
  kap<-if(pe<1)(po-pe)/(1-pe) else NA
  c(Acc=acc,BalAcc=if(!is.na(se)&&!is.na(sp))(se+sp)/2 else NA,Sens=se,Spec=sp,F1=f1,MCC=mcc,Kappa=kap)}
.youden <- function(y,p){if(length(unique(y))<2L)return(0.5)
  u<-sort(unique(p));if(length(u)==1L)return(0.5)
  eps<-max(1e-12,sqrt(.Machine$double.eps))
  cand<-c(min(u)-eps,(u[-1L]+u[-length(u)])/2,max(u)+eps)
  jstat<-vapply(cand,function(t){m<-.thr_metrics(y,p,t)
    if(is.na(m["Sens"])||is.na(m["Spec"]))-Inf else m["Sens"]+m["Spec"]-1},0)
  opt<-cand[jstat==max(jstat)]
  opt[which.min(abs(opt-0.5))]}                             
.platt_fit <- function(y,p){
  lo<-stats::qlogis(pmin(pmax(p,1e-6),1-1e-6))
  f<-tryCatch(stats::glm(y~lo,family=stats::binomial()),error=function(e)NULL)
  if(is.null(f)||!isTRUE(f$converged)||any(!is.finite(stats::coef(f))))stop("Platt calibration failed")
  f}
.platt_predict <- function(f,p){
  lo<-stats::qlogis(pmin(pmax(p,1e-6),1-1e-6))
  out<-as.vector(stats::predict(f,newdata=data.frame(lo=lo),type="response"))
  if(length(out)!=length(p)||any(!is.finite(out)))stop("Platt prediction failed");out}
.platt_oof <- function(y,p_oof,K=3){
  fo<-.strat_folds(y,K);out<-numeric(length(y))
  for(f in 1:K){tri<-which(fo!=f);vai<-which(fo==f)
    if(length(unique(y[tri]))<2){out[vai]<-mean(y[tri]);next}
    cf<-.platt_fit(y[tri],p_oof[tri]);out[vai]<-.platt_predict(cf,p_oof[vai])}
  out}

evaluate <- function(y_te,p_te,thr_oof){
  cs<-.cal(y_te,p_te);m50<-.thr_metrics(y_te,p_te,0.5);mY<-.thr_metrics(y_te,p_te,thr_oof)
  out<-c(Test_AUC=.auc(y_te,p_te),Test_PRAUC=.pr_auc(y_te,p_te),
    Test_Brier=mean((p_te-y_te)^2),Test_LogLoss=.ll(y_te,p_te),Test_ECE=.ece(y_te,p_te),
    Test_ECEadapt=.ece_adaptive(y_te,p_te),Test_ICI=.ici(y_te,p_te),
    Test_CalSlope=cs[1],Test_CalIntercept=cs[2],Thr_OOF=thr_oof,
    stats::setNames(m50,paste0("Test_",names(m50),"_50")),
    stats::setNames(mY,paste0("Test_",names(mY),"_Y")))
  names(out)<-sub("^Test_CalSlope.*","Test_CalSlope",names(out))
  names(out)<-sub("^Test_CalIntercept.*","Test_CalIntercept",names(out));out}
.thr_metrics_yhat <- function(y,yh){
  tp<-sum(yh==1&y==1);tn<-sum(yh==0&y==0);fp<-sum(yh==1&y==0);fn<-sum(yh==0&y==1)
  se<-if(tp+fn>0)tp/(tp+fn) else NA;sp<-if(tn+fp>0)tn/(tn+fp) else NA
  pr<-if(tp+fp>0)tp/(tp+fp) else NA
  f1<-if(!is.na(pr)&&!is.na(se)&&(pr+se)>0)2*pr*se/(pr+se) else NA
  den<-sqrt(as.numeric(tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  mcc<-if(den>0)(as.numeric(tp)*tn-as.numeric(fp)*fn)/den else NA
  acc<-(tp+tn)/length(y);pe<-((tp+fp)*(tp+fn)+(tn+fn)*(tn+fp))/length(y)^2
  kap<-if(pe<1)(acc-pe)/(1-pe) else NA
  c(Acc=acc,BalAcc=if(!is.na(se)&&!is.na(sp))(se+sp)/2 else NA,Sens=se,Spec=sp,F1=f1,MCC=mcc,Kappa=kap)}
evaluate_pooled <- function(y,p,thr_vec){
  cs<-.cal(y,p);m50<-.thr_metrics(y,p,0.5);mY<-.thr_metrics_yhat(y,as.integer(p>=thr_vec))
  out<-c(Test_AUC=.auc(y,p),Test_PRAUC=.pr_auc(y,p),
    Test_Brier=mean((p-y)^2),Test_LogLoss=.ll(y,p),Test_ECE=.ece(y,p),
    Test_ECEadapt=.ece_adaptive(y,p),Test_ICI=.ici(y,p),
    Test_CalSlope=cs[1],Test_CalIntercept=cs[2],Thr_OOF=stats::median(thr_vec),
    stats::setNames(m50,paste0("Test_",names(m50),"_50")),
    stats::setNames(mY,paste0("Test_",names(mY),"_Y")))
  names(out)<-sub("^Test_CalSlope.*","Test_CalSlope",names(out))
  names(out)<-sub("^Test_CalIntercept.*","Test_CalIntercept",names(out));out}
METRIC_NAMES <- c("Test_AUC","Test_PRAUC","Test_Brier","Test_LogLoss","Test_ECE",
  "Test_ECEadapt","Test_ICI","Test_CalSlope","Test_CalIntercept","Thr_OOF",
  paste0("Test_",c("Acc","BalAcc","Sens","Spec","F1","MCC","Kappa"),"_50"),
  paste0("Test_",c("Acc","BalAcc","Sens","Spec","F1","MCC","Kappa"),"_Y"))

.df <- function(X,cat_mask){d<-as.data.frame(X);names(d)<-paste0("V",seq_len(ncol(X)))
  for(j in which(cat_mask))d[[j]]<-factor(d[[j]]);d}
.dummy <- function(Xtr,Xte,cat_mask){
  a<-.df(Xtr,cat_mask);b<-.df(Xte,cat_mask)
  for(j in which(cat_mask)){lv<-levels(a[[j]]);b[[j]]<-factor(as.character(b[[j]]),levels=lv)}
  mm<-stats::model.matrix(~.-1,data=a);mt<-stats::model.matrix(~.-1,data=b)
  stopifnot(ncol(mm)==ncol(mt))
  cn<-sprintf("x%02d",seq_len(ncol(mm)))
  colnames(mm)<-cn;colnames(mt)<-cn
  list(tr=mm,te=mt)}


.attach_brcm_diag <- function(p,f){
  attr(p,"fit_diag")<-list(converged=as.integer(f$converged),iter=f$iter,lambda=f$lambda,
    n_inter=f$n_inter,n_screened_pairs=if(is.null(f$n_screened_pairs))NA_integer_ else f$n_screened_pairs,
    interaction_search_enabled=if(is.null(f$interaction_search_enabled))NA else f$interaction_search_enabled); p}
m_brcm <- function(Xtr,ytr,Xte,cm){f<-brcm_fit(Xtr,ytr,cm)
  if(!isTRUE(f$converged))stop("B-RCM++ did not converge")
  p<-brcm_predict(f,Xte);.attach_brcm_diag(p,f)}
m_abl_lin <- function(Xtr,ytr,Xte,cm){f<-brcm_fit(Xtr,ytr,cm,use_rank=FALSE,use_bernstein=FALSE,use_spearman=FALSE,use_inter=FALSE);if(!isTRUE(f$converged))stop("ablation did not converge");.attach_brcm_diag(brcm_predict(f,Xte),f)}
m_abl_rank <- function(Xtr,ytr,Xte,cm){f<-brcm_fit(Xtr,ytr,cm,use_bernstein=FALSE,use_spearman=FALSE,use_inter=FALSE);if(!isTRUE(f$converged))stop("ablation did not converge");.attach_brcm_diag(brcm_predict(f,Xte),f)}
m_abl_bern <- function(Xtr,ytr,Xte,cm){f<-brcm_fit(Xtr,ytr,cm,use_spearman=FALSE,use_inter=FALSE);if(!isTRUE(f$converged))stop("ablation did not converge");.attach_brcm_diag(brcm_predict(f,Xte),f)}
m_abl_spw <- function(Xtr,ytr,Xte,cm){f<-brcm_fit(Xtr,ytr,cm,use_inter=FALSE);if(!isTRUE(f$converged))stop("ablation did not converge");.attach_brcm_diag(brcm_predict(f,Xte),f)}
m_lr <- function(Xtr,ytr,Xte,cm){a<-.df(Xtr,cm);b<-.df(Xte,cm)
  for(j in which(cm))b[[j]]<-factor(as.character(b[[j]]),levels=levels(a[[j]]))
  f<-suppressWarnings(stats::glm(ytr~.,data=a,family=stats::binomial()))
  if(!isTRUE(f$converged))stop("LR (glm) did not converge")
  as.vector(suppressWarnings(stats::predict(f,newdata=b,type="response")))}
m_plr <- function(Xtr,ytr,Xte,cm){dd<-.dummy(Xtr,Xte,cm)
  fo<-.strat_folds(ytr,3);grid<-c(0.25,0.5,1,2,4,8,16)
  cv<-vapply(grid,function(lam){e<-c();ok<-TRUE
    for(f in 1:3){tri<-which(fo!=f);vai<-which(fo==f);if(length(unique(ytr[tri]))<2){ok<-FALSE;break}
      mu_f<-colMeans(dd$tr[tri,,drop=FALSE]);sd_f<-apply(dd$tr[tri,,drop=FALSE],2,stats::sd);sd_f[sd_f<1e-12]<-1
      Ztri<-scale(dd$tr[tri,,drop=FALSE],center=mu_f,scale=sd_f)      
      Zval<-scale(dd$tr[vai,,drop=FALSE],center=mu_f,scale=sd_f)
      ft<-.irls(cbind(1,Ztri),ytr[tri],c(0,rep(lam,ncol(Ztri))))
      if(!isTRUE(ft$converged)){ok<-FALSE;break}                     
      e<-c(e,.ll(ytr[vai],1/(1+exp(-drop(cbind(1,Zval)%*%ft$beta)))))}
    if(ok&&length(e)==3L)mean(e) else NA_real_},numeric(1))          
  lam<-if(all(is.na(cv)))1 else grid[which.min(cv)]
  mu<-colMeans(dd$tr);sdv<-apply(dd$tr,2,stats::sd);sdv[sdv<1e-12]<-1  
  Ztr<-scale(dd$tr,center=mu,scale=sdv);Zte<-scale(dd$te,center=mu,scale=sdv)
  ft<-.irls(cbind(1,Ztr),ytr,c(0,rep(lam,ncol(Ztr))))
  if(!isTRUE(ft$converged))stop("PLR IRLS did not converge")
  as.vector(1/(1+exp(-drop(cbind(1,Zte)%*%ft$beta))))}
m_gam <- function(Xtr,ytr,Xte,cm){a<-.df(Xtr,cm);b<-.df(Xte,cm)
  for(j in which(cm))b[[j]]<-factor(as.character(b[[j]]),levels=levels(a[[j]]))
  cont<-which(!cm);nu<-vapply(cont,function(j)length(unique(Xtr[,j])),integer(1))
  smooth_idx<-cont[nu>=5];linear_idx<-cont[nu>=2 & nu<5];cat_idx<-which(cm)  
  vname<-function(idx)if(length(idx))paste0("V",idx) else character(0)       
  smooth_terms<-if(length(smooth_idx))vapply(smooth_idx,function(j)sprintf("s(V%d,k=%d)",j,min(5L,length(unique(Xtr[,j]))-1L)),character(1)) else character(0)
  terms<-c(smooth_terms,vname(linear_idx),vname(cat_idx))
  fo<-if(length(terms))stats::as.formula(paste("ytr ~",paste(terms,collapse="+"))) else stats::as.formula("ytr ~ 1")
  f<-tryCatch(mgcv::gam(fo,data=a,family=stats::binomial(),method="REML"),error=function(e)NULL)
  if(is.null(f))stop("GAM fit failed")
  p<-as.vector(mgcv::predict.gam(f,newdata=b,type="response"))
  if(length(p)!=nrow(Xte)||any(!is.finite(p)))stop("GAM produced invalid predictions");p}
m_dt <- function(Xtr,ytr,Xte,cm){a<-.df(Xtr,cm);b<-.df(Xte,cm)
  for(j in which(cm))b[[j]]<-factor(as.character(b[[j]]),levels=levels(a[[j]]))
  f<-rpart::rpart(factor(ytr)~.,data=a,method="class",control=rpart::rpart.control(cp=0.001,xval=10))
  cp<-f$cptable;best<-cp[which.min(cp[,"xerror"]),"CP"]
  f<-rpart::prune(f,cp=best);as.vector(stats::predict(f,newdata=b,type="prob")[,2])}
m_rf <- function(Xtr,ytr,Xte,cm){a<-.df(Xtr,cm);b<-.df(Xte,cm)
  for(j in which(cm))b[[j]]<-factor(as.character(b[[j]]),levels=levels(a[[j]]))
  f<-randomForest::randomForest(x=a,y=factor(ytr),ntree=500,mtry=max(1,floor(sqrt(ncol(a)))))
  as.vector(stats::predict(f,newdata=b,type="prob")[,2])}
m_svm <- function(Xtr,ytr,Xte,cm){dd<-.dummy(Xtr,Xte,cm)
  f<-e1071::svm(x=dd$tr,y=factor(ytr),kernel="radial",cost=1,probability=TRUE,scale=TRUE)
  pr<-attr(stats::predict(f,dd$te,probability=TRUE),"probabilities")
  as.vector(pr[,which(colnames(pr)=="1")])}
m_xgb <- function(Xtr,ytr,Xte,cm){dd<-.dummy(Xtr,Xte,cm)
  f<-xgboost::xgboost(data=dd$tr,label=ytr,nrounds=100,eta=0.1,max_depth=3,
    objective="binary:logistic",verbose=0,nthread=1)
  as.vector(stats::predict(f,dd$te))}
m_lgbm <- function(Xtr,ytr,Xte,cm){dd<-.dummy(Xtr,Xte,cm)
  ds<-lightgbm::lgb.Dataset(data=dd$tr,label=ytr)
  f<-lightgbm::lgb.train(params=list(objective="binary",learning_rate=0.1,num_leaves=15,
    verbosity=-1L,num_threads=1L),data=ds,nrounds=100)
  as.vector(stats::predict(f,dd$te))}
m_catb <- function(Xtr,ytr,Xte,cm){a<-.df(Xtr,cm);b<-.df(Xte,cm)
  for(j in which(cm))b[[j]]<-factor(as.character(b[[j]]),levels=levels(a[[j]]))
  catf<-which(cm)-1L                                   
  pl<-catboost::catboost.load_pool(data=a,label=ytr,cat_features=catf)
  pt<-catboost::catboost.load_pool(data=b,cat_features=catf)
  f<-catboost::catboost.train(pl,params=list(loss_function="Logloss",iterations=100,
    learning_rate=0.1,depth=6,logging_level="Silent",thread_count=1))
  as.vector(catboost::catboost.predict(f,pt,prediction_type="Probability"))}

.ebm_bind_python <- function(){
  if(!requireNamespace("reticulate", quietly=TRUE)) return("reticulate not installed")
  if(isTRUE(tryCatch(reticulate::py_available(initialize=FALSE), error=function(e) FALSE)))
    return("already bound")                      
  want <- Sys.getenv("BRCM_PYTHON", Sys.getenv("RETICULATE_PYTHON", ""))
  if(nzchar(want) && file.exists(want)){
    ok <- tryCatch({reticulate::use_python(want, required=TRUE); TRUE}, error=function(e) FALSE)
    if(ok) return(paste("BRCM_PYTHON:", want))}
  envs <- tryCatch(reticulate::conda_list(), error=function(e) NULL)
  if(!is.null(envs) && nrow(envs)){
    hit <- envs$python[envs$name == "brcm-ebm"]
    if(length(hit) && file.exists(hit[1])){
      ok <- tryCatch({reticulate::use_condaenv("brcm-ebm", required=TRUE); TRUE}, error=function(e) FALSE)
      if(ok) return("conda env: brcm-ebm")}
    for(k in seq_len(nrow(envs))){
      pp <- envs$python[k]; if(!file.exists(pp)) next
      lib <- file.path(dirname(pp), "Lib", "site-packages", "interpret")
      lib2 <- file.path(dirname(dirname(pp)), "lib", "site-packages", "interpret")
      if(dir.exists(lib) || dir.exists(lib2)){
        ok <- tryCatch({reticulate::use_python(pp, required=TRUE); TRUE}, error=function(e) FALSE)
        if(ok) return(paste("conda env:", envs$name[k]))}}}
  "reticulate default"}

.EBM_BIND <- tryCatch(.ebm_bind_python(), error=function(e) paste("bind failed:", conditionMessage(e)))

.EBM_OK <- local({
  if(!requireNamespace("reticulate", quietly=TRUE)) return(FALSE)
  if(!isTRUE(tryCatch(reticulate::py_module_available("interpret"),
                      error=function(e) FALSE))) return(FALSE)
  res <- tryCatch({
    gl <- reticulate::import("interpret.glassbox", delay_load=FALSE)
    np <- reticulate::import("numpy", delay_load=FALSE)
    set.seed(1); Xc <- matrix(stats::rnorm(120), 60, 2)
    yc <- as.integer(Xc[,1] + stats::rnorm(60) > 0)
    if(length(unique(yc)) < 2L) yc[1] <- 1L - yc[1]
    cc <- gl$ExplainableBoostingClassifier(interactions=1L, random_state=1L, n_jobs=1L)
    cc$fit(np$asarray(Xc), np$asarray(yc))
    pp <- reticulate::py_to_r(cc$predict_proba(np$asarray(Xc)))
    is.matrix(pp) && ncol(pp) == 2L && all(is.finite(pp))
  }, error=function(e){assign(".EBM_ERR", conditionMessage(e), envir=.GlobalEnv); FALSE})
  isTRUE(res) })

brcm_ebm_check <- function(){
  cat("Python binding : ", .EBM_BIND, "\n", sep="")
  if(requireNamespace("reticulate", quietly=TRUE)){
    cfg <- tryCatch(reticulate::py_config(), error=function(e) NULL)
    if(!is.null(cfg)) cat("Interpreter    : ", cfg$python, "\n", sep="")
    cat("interpret found: ",
        isTRUE(tryCatch(reticulate::py_module_available("interpret"), error=function(e) FALSE)),
        "\n", sep="")}
  cat("Trial fit      : ", .EBM_OK, "\n", sep="")
  if(!.EBM_OK && exists(".EBM_ERR", envir=.GlobalEnv))
    cat("Error          : ", get(".EBM_ERR", envir=.GlobalEnv), "\n", sep="")
  invisible(.EBM_OK)}

m_ebm <- function(Xtr,ytr,Xte,cm){
  if(!.EBM_OK) stop("EBM unavailable: reticulate or the Python interpret package is not installed")
  dd <- .dummy(Xtr,Xte,cm)                     
  gl <- reticulate::import("interpret.glassbox", delay_load=FALSE)
  np <- reticulate::import("numpy", delay_load=FALSE)
  clf <- gl$ExplainableBoostingClassifier(interactions=4L, random_state=1L, n_jobs=1L)
  clf$fit(np$asarray(dd$tr), np$asarray(as.integer(ytr)))
  pr <- clf$predict_proba(np$asarray(dd$te))
  as.vector(reticulate::py_to_r(pr)[,2])}



build_registry <- function(){
  reg <- list("B-RCM++"=m_brcm,"LR"=m_lr,"PLR"=m_plr,"DT"=m_dt)
  if(has("mgcv")) reg[["GAM"]] <- m_gam
  if(has("randomForest")) reg[["RF"]] <- m_rf
  if(has("e1071")) reg[["SVM"]] <- m_svm
  if(has("xgboost")) reg[["XGBoost"]] <- m_xgb
  if(has("lightgbm")) reg[["LightGBM"]] <- m_lgbm
  if(has("catboost")) reg[["CatBoost"]] <- m_catb
  if(.EBM_OK){ reg[["EBM"]] <- m_ebm
  } else if(isTRUE(STRICT_PKGS)) stop(
    "EBM is required for the reported run but the Python 'interpret' stack is not usable.\n",
    "  Install it with:  reticulate::py_install(\"interpret\", pip = TRUE)\n",
    "  Then restart the R session, because reticulate binds to an interpreter once per session.\n",
    "  If the module imports but the fit fails with a TerminatedWorkerError, the cause is\n",
    "  joblib competing with the R cluster; this file already passes n_jobs = 1 to avoid it.\n",
    "  Run brcm_ebm_check() to see which interpreter was bound and what the trial fit reported.\n",
    "  Set BRCM_STRICT=FALSE for an exploratory run without EBM.")
  reg[["A1-LinearRidge"]]<-m_abl_lin; reg[["A2-RankLinear"]]<-m_abl_rank
  reg[["A3-RankBernstein"]]<-m_abl_bern; reg[["A4-PlusSpearman"]]<-m_abl_spw
  reg }
CALIBRATE <- c("RF","XGBoost","LightGBM","CatBoost")  
oof_probs_rd <- function(fitter,Xraw,ytr,cm,K=OOF_K){
  fo<-.strat_folds(ytr,K);p<-numeric(length(ytr))
  for(f in 1:K){itr<-which(fo!=f);iva<-which(fo==f)
    if(length(unique(ytr[itr]))<2){p[iva]<-mean(ytr[itr]);next}
    fillf<-.impute_fit(Xraw[itr,,drop=FALSE],cm)             
    Xi<-.impute_apply(Xraw[itr,,drop=FALSE],fillf)
    Xv<-.impute_apply(Xraw[iva,,drop=FALSE],fillf)
    pv<-tryCatch(fitter(Xi,ytr[itr],Xv,cm),error=function(e)e)
    if(inherits(pv,"error")||length(pv)!=length(iva)||any(!is.finite(pv)))stop("RD OOF inner fit failed")
    p[iva]<-pv}
  p}

.as_num <- function(v){if(is.factor(v))as.numeric(v) else if(is.character(v))as.numeric(factor(v)) else as.numeric(v)}
load_datasets <- function(){
  ds<-list();e<-new.env()
  utils::data(list=c("PimaIndiansDiabetes","BreastCancer","Ionosphere","Sonar","HouseVotes84"),
              package="mlbench",envir=e)
  d<-e$PimaIndiansDiabetes;X<-as.matrix(d[,1:8])
  for(v in c("glucose","pressure","triceps","insulin","mass")) X[X[,v]==0,v]<-NA
  ds[["PimaIndiansDiabetes"]]<-list(X=X,y=as.integer(d$diabetes=="pos"),
    cat=rep(FALSE,8),label="Pima Indians Diabetes",pos="pos",
    ptype="continuous (clinical)",miss="physiologically impossible zeros in glucose, pressure, triceps, insulin, mass recoded as missing, then median-imputed within training folds",enc="diabetes = pos")
  d<-e$BreastCancer;Xd<-d[,2:10];X<-vapply(Xd,.as_num,numeric(nrow(d)))
  ds[["BreastCancer"]]<-list(X=X,y=as.integer(d$Class=="malignant"),
    cat=rep(FALSE,ncol(X)),label="Breast Cancer (Wisconsin)",pos="malignant",
    ptype="ordinal 1-10",miss="native NAs median-imputed within training folds; Id column dropped",enc="class = malignant")
  d<-e$Ionosphere;Xd<-d[,1:34];keep<-vapply(Xd,function(v)length(unique(v[!is.na(v)]))>1,TRUE)
  Xd<-Xd[,keep,drop=FALSE];X<-vapply(Xd,.as_num,numeric(nrow(d)))
  ds[["Ionosphere"]]<-list(X=X,y=as.integer(d$Class=="bad"),
    cat=vapply(Xd,function(v)is.factor(v)||length(unique(v[!is.na(v)]))<=2,TRUE),
    label="Ionosphere",pos="bad",
    ptype="continuous (correlated)",miss="none; constant predictor V2 dropped",enc="class = bad")
  d<-e$Sonar;X<-as.matrix(d[,1:60])
  ds[["Sonar"]]<-list(X=X,y=as.integer(d$Class=="M"),cat=rep(FALSE,60),label="Sonar",pos="M",
    ptype="continuous (p large vs n)",miss="none",enc="class = M (mine)")
  d<-e$HouseVotes84;Xd<-d[,2:17];X<-vapply(Xd,.as_num,numeric(nrow(d)))
  ds[["HouseVotes84"]]<-list(X=X,y=as.integer(d$Class=="democrat"),
    cat=rep(TRUE,ncol(X)),label="Congressional Voting Records",pos="democrat",
    ptype="binary categorical",miss="abstentions/NAs mode-imputed within training folds",enc="party = democrat")
  ds}

.impute_fit <- function(X,cat_mask){
  vapply(seq_len(ncol(X)),function(j){v<-X[,j];v<-v[!is.na(v)]
    if(!length(v))return(0)
    if(cat_mask[j]){tb<-table(v);as.numeric(names(tb)[which.max(tb)])} else stats::median(v)},numeric(1))}
.impute_apply <- function(X,fill){for(j in seq_len(ncol(X))){i<-is.na(X[,j]);if(any(i))X[i,j]<-fill[j]};X}

.no_default_device <- function(){
  if(is.null(getOption("bitmapType")) || TRUE) options(device = function(...) grDevices::png(tempfile(fileext=".png")))
  invisible(NULL)}
.no_default_device()

.fmt_hms <- function(secs){
  if(!is.finite(secs) || secs < 0) return("--:--:--")
  secs <- as.numeric(secs)
  h <- floor(secs/3600); m <- floor((secs - 3600*h)/60); s <- floor(secs - 3600*h - 60*m)
  sprintf("%02d:%02d:%02d", h, m, s)}

.prog_new <- function(total, label = "", quiet = FALSE){
  pr <- list(total = as.integer(total), i = 0L, t0 = Sys.time(),
             label = label, quiet = isTRUE(quiet))
  if(!pr$quiet){
    cat(sprintf("\n[%s] %s: starting %d iterations\n",
                format(pr$t0, "%H:%M:%S"), label, pr$total))
    flush.console()}
  pr}

.prog_tick <- function(pr, note = ""){
  pr$i <- pr$i + 1L
  if(pr$quiet) return(pr)
  el  <- as.numeric(difftime(Sys.time(), pr$t0, units = "secs"))
  per <- if(pr$i > 0L) el/pr$i else NA_real_
  eta <- if(is.finite(per)) per * (pr$total - pr$i) else NA_real_
  cat(sprintf("[%s] %s %d/%d (%5.1f%%) | elapsed %s | eta %s%s\n",
              format(Sys.time(), "%H:%M:%S"), pr$label, pr$i, pr$total,
              100 * pr$i / pr$total, .fmt_hms(el), .fmt_hms(eta),
              if(nzchar(note)) paste0(" | ", note) else ""))
  flush.console()
  pr}

.prog_done <- function(pr){
  if(pr$quiet) return(invisible(NULL))
  el <- as.numeric(difftime(Sys.time(), pr$t0, units = "secs"))
  cat(sprintf("[%s] %s: finished %d/%d in %s\n",
              format(Sys.time(), "%H:%M:%S"), pr$label, pr$i, pr$total, .fmt_hms(el)))
  flush.console()
  invisible(NULL)}

.cv_folds <- function(y,K,seed){                                   
  old<-if(exists(".Random.seed",envir=.GlobalEnv))get(".Random.seed",envir=.GlobalEnv) else NULL
  set.seed(seed);f<-integer(length(y))
  for(cl in unique(y)){i<-which(y==cl);f[i]<-sample(rep(seq_len(K),length.out=length(i)))}
  if(!is.null(old))assign(".Random.seed",old,envir=.GlobalEnv);f}

run_dataset <- function(ds,registry,n_repeats=N_REPEATS,n_folds=N_FOLDS,base_seed=BASE_SEED,cl=NULL,
                        progress=TRUE,progress_label=""){
  X0<-ds$X;y<-ds$y;cm<-ds$cat;n<-length(y)
  methods<-names(registry)
  one_repeat<-function(rp){
    folds<-.cv_folds(y,n_folds,base_seed+rp)
    P<-list();TH<-list()
    for(nm in methods){for(v in c(nm,paste0(nm,"+Platt"))){P[[v]]<-rep(NA_real_,n);TH[[v]]<-rep(NA_real_,n)}}
    itr<-list()                                 
    for(fl in seq_len(n_folds)){
      tri<-which(folds!=fl);tei<-which(folds==fl)
      if(length(unique(y[tri]))<2)next
      Xtr_raw<-X0[tri,,drop=FALSE];Xte_raw<-X0[tei,,drop=FALSE];ytr<-y[tri]
      fill<-.impute_fit(Xtr_raw,cm)                               
      Xtr<-.impute_apply(Xtr_raw,fill);Xte<-.impute_apply(Xte_raw,fill)
      for(nm in methods){
        set.seed(base_seed*7L+rp*1000L+fl*50L+.name_seed(nm))      
        fitter<-registry[[nm]]
        r<-tryCatch({
          p_inner<-oof_probs_rd(fitter,Xtr_raw,ytr,cm)             
          thr<-.youden(ytr,p_inner)
          p_te<-fitter(Xtr,ytr,Xte,cm);if(any(!is.finite(p_te)))stop("non-finite")
          list(te=p_te,thr=thr,p_inner=p_inner)},error=function(e)NULL)
        if(is.null(r))next
        P[[nm]][tei]<-r$te;TH[[nm]][tei]<-r$thr
        if(nm=="B-RCM++"){fd<-attr(r$te,"fit_diag")
          if(!is.null(fd)) itr[[length(itr)+1L]]<-data.frame(repeat_id=rp,fold=fl,
            n_inter=as.integer(fd$n_inter),
            n_screened_pairs=as.integer(fd$n_screened_pairs),
            search_enabled=as.logical(fd$interaction_search_enabled),stringsAsFactors=FALSE)}
        if(nm %in% CALIBRATE){
          cal<-tryCatch({cf<-.platt_fit(ytr,r$p_inner)               
            thc<-.youden(ytr,.platt_oof(ytr,r$p_inner))              
            list(cal=.platt_predict(cf,r$te),thr=thc)},error=function(e)NULL)
          if(!is.null(cal)){P[[paste0(nm,"+Platt")]][tei]<-cal$cal;TH[[paste0(nm,"+Platt")]][tei]<-cal$thr}}}}
    out<-list()
    for(nm in names(P)){p<-P[[nm]];th<-TH[[nm]];ok<-is.finite(p)&is.finite(th)
      if(sum(ok)!=n||length(unique(y[ok]))<2L)next     
      out[[nm]]<-evaluate_pooled(y[ok],p[ok],th[ok])}
    list(id=paste0("r",rp),metrics=out,
         inter=if(length(itr))do.call(rbind,itr) else NULL)}
  if(isTRUE(progress) && n_repeats >= 2L && PROG_CHUNKS >= 1L){
    lbl<-if(nzchar(progress_label)) progress_label else sprintf("%s repeats",ds$label)
    t0<-Sys.time(); res<-vector("list",n_repeats)
    .newcl<-function(){
      k<-tryCatch(parallel::makeCluster(n_cores),error=function(e)NULL)
      if(is.null(k)) return(NULL)
      ok<-tryCatch({parallel::clusterExport(k,varlist=ls(envir=.GlobalEnv,all.names=TRUE),
                                            envir=.GlobalEnv)
                    parallel::clusterEvalQ(k,{suppressMessages({library(stats);library(utils)})})
                    if(isTRUE(.EBM_OK)) parallel::clusterEvalQ(k,{
                      try(.ebm_bind_python(), silent=TRUE)
                      try(reticulate::import("interpret.glassbox", delay_load=FALSE), silent=TRUE)
                      NULL})
                    TRUE},error=function(e)FALSE)
      if(!ok){try(parallel::stopCluster(k),silent=TRUE);return(NULL)}
      k}
    nw<-max(1L, if(is.null(cl)) 1L else length(cl))
    starts<-seq(1L, n_repeats, by=nw)
    for(st in starts){
      idx<-st:min(st+nw-1L, n_repeats)
      part<-tryCatch(if(!is.null(cl)) parallel::parLapply(cl,idx,one_repeat) else lapply(idx,one_repeat),
                     error=function(e)NULL)
      if(is.null(part)){
        cat(sprintf("\n    [!] repeats %d-%d of %s lost a worker; rebuilding and isolating\n",
                    min(idx), max(idx),
                    if(nzchar(progress_label)) progress_label else ds$label));flush.console()
        if(!is.null(cl)) try(parallel::stopCluster(cl),silent=TRUE)
        cl<-.newcl()
        part<-lapply(idx,function(u)
          tryCatch(if(!is.null(cl)) parallel::parLapply(cl,u,one_repeat)[[1]] else one_repeat(u),
                   error=function(e){cat(sprintf("    [!] repeat %d could not be completed\n",u));NULL}))}
      res[idx]<-part
      el<-as.numeric(difftime(Sys.time(),t0,units="secs"))
      dn<-max(idx); eta<-el/dn*(n_repeats-dn)
      cat(sprintf("\r    %s %d/%d (%5.1f%%) | elapsed %s | eta %s   ",
                  lbl, dn, n_repeats, 100*dn/n_repeats, .fmt_hms(el), .fmt_hms(eta)))
      flush.console()}
    cat("\n")
  } else {
    res<-tryCatch(if(!is.null(cl))parallel::parLapply(cl,seq_len(n_repeats),one_repeat) else lapply(seq_len(n_repeats),one_repeat),
                  error=function(e){cat("    [!] dispatch failed:",conditionMessage(e),"\n");vector("list",n_repeats)})
  }
  allm<-unique(unlist(lapply(res,function(z)names(z$metrics))))
  raw<-list()
  for(nm in allm){ids<-character(0);rows<-list()
    for(z in res)if(!is.null(z$metrics[[nm]])){ids<-c(ids,z$id);rows[[length(rows)+1]]<-z$metrics[[nm]]}
    if(!length(rows))next
    M<-do.call(rbind,rows);rownames(M)<-ids                       
    raw[[nm]]<-M}
  diag_methods<-unique(c(methods,paste0(intersect(methods,CALIBRATE),"+Platt")))
  dg<-do.call(rbind,lapply(diag_methods,function(nm){
    att<-n_repeats;okn<-if(!is.null(raw[[nm]]))nrow(raw[[nm]]) else 0L
    data.frame(Method=nm,Attempted=att,Successful=okn,
      Failure_rate=round(100*(att-okn)/att,2),stringsAsFactors=FALSE)}))
  keep<-!vapply(res,is.null,TRUE)
  if(!all(keep)){cat(sprintf("    [!] %d of %d repeats lost on this dataset\n",sum(!keep),n_repeats))
    res<-res[keep]}
  itab<-do.call(rbind,lapply(res,function(z)z$inter))
  list(raw=raw,diag=dg,inter=itab,cl=cl)}

.rank_biserial <- function(d){d<-d[is.finite(d)];d<-d[d!=0];if(length(d)<2)return(NA_real_)
  r<-rank(abs(d));(sum(r[d>0])-sum(r[d<0]))/sum(r)}
LOWER_BETTER <- c(Test_AUC=FALSE,Test_PRAUC=FALSE,Test_Brier=TRUE,Test_LogLoss=TRUE,
  Test_ECE=TRUE,Test_ECEadapt=TRUE,Test_ICI=TRUE)
paired_analysis <- function(raw_all,reference="B-RCM++",metrics=names(LOWER_BETTER)){
  rows<-list();k<-0
  for(dsn in names(raw_all)){R<-raw_all[[dsn]];ref<-R[[reference]];if(is.null(ref))next
    for(cmp in setdiff(names(R),reference)){oth<-R[[cmp]];if(is.null(oth))next
      ids<-intersect(rownames(ref),rownames(oth));if(length(ids)<3)next
      for(met in metrics){
        if(!(met%in%colnames(ref))||!(met%in%colnames(oth)))next
        d<-ref[ids,met]-oth[ids,met];d<-d[is.finite(d)];if(length(d)<3)next
        dd<-if(LOWER_BETTER[[met]]) -d else d
        wt<-tryCatch(stats::wilcox.test(dd,conf.int=TRUE,exact=FALSE),error=function(e)NULL)
        k<-k+1
        rows[[k]]<-data.frame(Dataset=dsn,Comparator=cmp,Metric=met,n_common=length(ids),
          HL=if(!is.null(wt))unname(wt$estimate) else stats::median(dd),
          CI_lo=if(!is.null(wt))wt$conf.int[1] else NA_real_,
          CI_hi=if(!is.null(wt))wt$conf.int[2] else NA_real_,
          rank_biserial=.rank_biserial(dd),median_diff=stats::median(dd),
          wins=sum(dd>0),ties=sum(dd==0),losses=sum(dd<0),
          p_raw=if(!is.null(wt))wt$p.value else NA_real_,stringsAsFactors=FALSE)}}}
  D<-do.call(rbind,rows);if(is.null(D))return(NULL)
  D$p_BH<-NA_real_
  for(cmp in unique(D$Comparator))for(met in unique(D$Metric)){
    i<-which(D$Comparator==cmp&D$Metric==met);if(length(i))D$p_BH[i]<-stats::p.adjust(D$p_raw[i],method="BH")}
  D}

.agg <- function(raw_all){rows<-list();k<-0
  for(dsn in names(raw_all))for(nm in names(raw_all[[dsn]])){M<-raw_all[[dsn]][[nm]]
    for(met in colnames(M)){v<-M[,met];v<-v[is.finite(v)];if(!length(v))next
      k<-k+1;rows[[k]]<-data.frame(Dataset=dsn,Model=nm,Metric=met,N=length(v),
        Mean=mean(v),MC_SE=stats::sd(v)/sqrt(length(v)),SD=stats::sd(v),Median=stats::median(v),
        Q1=stats::quantile(v,.25,names=FALSE),Q3=stats::quantile(v,.75,names=FALSE),
        stringsAsFactors=FALSE)}}
  do.call(rbind,rows)}

write_tables_rd <- function(res,out_dir=OUT_DIR){
  td<-file.path(out_dir,"tables");dir.create(td,recursive=TRUE,showWarnings=FALSE)
  S<-.agg(res$raw);utils::write.csv(S,file.path(td,"RD_summary.csv"),row.names=FALSE)
  models<-unique(S$Model)
  fmt<-function(md,q1,q3,pct)if(is.na(md))"--" else if(pct)sprintf("%.1f (%.1f)",100*md,100*(q3-q1)) else sprintf("%.3f (%.3f)",md,q3-q1)
  METS<-list(c("Test_AUC","AUC","1"),c("Test_PRAUC","PR-AUC","1"),c("Test_Brier","Brier","0"),
    c("Test_LogLoss","Log-loss","0"),c("Test_ICI","ICI","0"),
    c("Test_ECEadapt","ECE (adaptive)","0"),c("Test_ECE","ECE (equal-width)","0"),
    c("Test_CalSlope","Cal. slope","0"),
    c("Test_BalAcc_Y","Bal. acc.","1"),c("Test_F1_Y","F1","1"),c("Test_MCC_Y","MCC","0"))
  ti<-1
  for(dsn in names(res$raw)){blk<-list()
    for(me in METS){r<-data.frame(Metric=me[2],stringsAsFactors=FALSE)
      for(m in models){x<-S[S$Dataset==dsn&S$Model==m&S$Metric==me[1],]
        r[[m]]<-if(nrow(x))fmt(x$Median[1],x$Q1[1],x$Q3[1],me[3]=="1") else "--"}
      blk[[length(blk)+1]]<-r}
    utils::write.csv(do.call(rbind,blk),file.path(td,sprintf("RT%d_%s.csv",ti,dsn)),row.names=FALSE);ti<-ti+1}
  o<-do.call(rbind,lapply(models,function(m){r<-S[S$Model==m,]
    g<-function(met)mean(r$Mean[r$Metric==met],na.rm=TRUE)
    data.frame(Model=m,AUC=g("Test_AUC"),PR_AUC=g("Test_PRAUC"),Brier=g("Test_Brier"),
      LogLoss=g("Test_LogLoss"),ICI=g("Test_ICI"),ECEadapt=g("Test_ECEadapt"),
      ECE=g("Test_ECE"),BalAcc_Y=g("Test_BalAcc_Y"),stringsAsFactors=FALSE)}))
  au<-S[S$Metric=="Test_AUC",];w<-stats::reshape(au[,c("Dataset","Model","Mean")],idvar="Dataset",
    timevar="Model",direction="wide");rk<-t(apply(-as.matrix(w[,-1]),1,rank,na.last="keep"))
  colnames(rk)<-sub("^Mean\\.","",colnames(w)[-1]);o$MeanRank_AUC<-round(colMeans(rk,na.rm=TRUE)[o$Model],2)
  utils::write.csv(o,file.path(td,"RT_overall.csv"),row.names=FALSE)
  if(!is.null(res$paired)){P<-res$paired
    P$cell<-sprintf("%+.4f%s (%+.4f, %+.4f) [%+.2f]",P$HL,ifelse(!is.na(P$p_BH)&P$p_BH<0.05,"*",""),
      P$CI_lo,P$CI_hi,P$rank_biserial)
    utils::write.csv(P,file.path(td,"RD_paired.csv"),row.names=FALSE)}
  utils::write.csv(res$diag,file.path(td,"RD_diagnostics.csv"),row.names=FALSE)
  if(!is.null(res$inter_summary)) utils::write.csv(res$inter_summary,file.path(td,"RD_interaction_summary.csv"),row.names=FALSE)
  if(!is.null(res$inter_fold))    utils::write.csv(res$inter_fold,file.path(td,"RD_interaction_by_fold.csv"),row.names=FALSE)
  if(has("openxlsx")){wb<-openxlsx::createWorkbook()
    ad<-function(n,d){openxlsx::addWorksheet(wb,n);openxlsx::writeData(wb,n,d)}
    ad("Overall",o);ad("Summary_long",S);if(!is.null(res$paired))ad("Paired",res$paired)
    ad("Diagnostics",res$diag)
    openxlsx::saveWorkbook(wb,file.path(out_dir,"RealData_Tables.xlsx"),overwrite=TRUE)}
  invisible(td)}

write_figures_rd <- function(res,datasets,out_dir=OUT_DIR){
  fd<-file.path(out_dir,"figures");dir.create(fd,recursive=TRUE,showWarnings=FALSE)
  S<-.agg(res$raw)
  core<-intersect(c("B-RCM++","EBM","LR","PLR","GAM","DT","RF","SVM","XGBoost","LightGBM","CatBoost"),unique(S$Model))
  cols<-grDevices::hcl.colors(length(core),"Dark 3");dsn<-names(res$raw)
  for(mm in c("Test_AUC","Test_LogLoss")){
    png(file.path(fd,paste0("RF_",sub("Test_","",mm),"_by_dataset.png")),width=3000,height=1800,res=300)
    M<-vapply(core,function(m)vapply(dsn,function(d)mean(S$Mean[S$Model==m&S$Dataset==d&S$Metric==mm],na.rm=TRUE),0),numeric(length(dsn)))
    op<-par(mar=c(9,4,3,1));barplot(t(M),beside=TRUE,col=cols,las=2,
      ylab=paste("Mean",sub("Test_","",mm)),main=paste(sub("Test_","",mm),"by dataset"))
    legend("topright",legend=core,fill=cols,cex=.6,bty="n");par(op);dev.off()}
  if(!is.null(res$paired)){P<-res$paired[res$paired$Metric=="Test_AUC",]
    if(nrow(P)){png(file.path(fd,"RF_HL_forest_AUC.png"),width=2400,height=2800,res=300)
      P<-P[order(P$HL),];op<-par(mar=c(5,10,3,2))
      plot(P$HL,seq_len(nrow(P)),pch=16,cex=.6,yaxt="n",ylab="",
        xlab="Hodges-Lehmann difference in AUC (positive favours B-RCM++)",
        main="Paired effects by dataset and competitor",xlim=range(c(P$CI_lo,P$CI_hi),na.rm=TRUE))
      segments(P$CI_lo,seq_len(nrow(P)),P$CI_hi,seq_len(nrow(P)),col="grey40")
      axis(2,at=seq_len(nrow(P)),labels=paste(P$Dataset,P$Comparator,sep=" | "),las=2,cex.axis=.45)
      abline(v=0,lty=2,col="red");par(op);dev.off()}}
  d<-datasets[["PimaIndiansDiabetes"]]
  fill<-.impute_fit(d$X,d$cat);Xi<-.impute_apply(d$X,fill)
  f<-tryCatch(brcm_fit(Xi,d$y,d$cat),error=function(e)NULL)
  if(!is.null(f)){brcm_plot_shapes(f,colnames(d$X),file=file.path(fd,"RF_pima_components.png"),ncol_grid=3)
    if(length(f$inter))brcm_plot_interaction(f,1,colnames(d$X),file=file.path(fd,"RF_pima_interaction.png"))}
  invisible(fd)}

run_full_realdata <- function(n_repeats=N_REPEATS,n_folds=N_FOLDS,n_cores=N_CORES,
                              out_dir=OUT_DIR,base_seed=BASE_SEED,which_ds=NULL){
  .check_pkgs();dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  datasets<-load_datasets();if(!is.null(which_ds))datasets<-datasets[which_ds]
  registry<-build_registry()
  cat("Datasets:",paste(names(datasets),collapse=", "),"\n")
  cat("Methods :",paste(names(registry),collapse=", "),"\n")
  n_cores<-max(1L,min(as.integer(n_cores),as.integer(n_repeats)))   
  cl<-NULL
  if(n_cores>1&&has("parallel")){cl<-parallel::makeCluster(n_cores)
    parallel::clusterExport(cl,varlist=ls(envir=.GlobalEnv,all.names=TRUE),envir=.GlobalEnv)
    parallel::clusterEvalQ(cl,{suppressMessages({
      for(p in c("mgcv","rpart","randomForest","e1071","xgboost","lightgbm","catboost"))
        if(requireNamespace(p,quietly=TRUE))library(p,character.only=TRUE)});NULL})
    if(isTRUE(.EBM_OK)){
      try(parallel::clusterEvalQ(cl,{
        try(.ebm_bind_python(), silent=TRUE)
        try(reticulate::import("interpret.glassbox", delay_load=FALSE), silent=TRUE)
        NULL}), silent=TRUE)
      chk<-tryCatch(unlist(parallel::clusterEvalQ(cl,
             isTRUE(tryCatch(reticulate::py_module_available("interpret"),
                             error=function(e) FALSE)))),
           error=function(e) rep(FALSE, length(cl)))
      if(!all(chk)){
        msg <- sprintf(
          "EBM is available in the master process but only %d of %d workers could reach it; EBM results would be incomplete.",
          sum(chk), length(chk))
        if(isTRUE(STRICT_PKGS)) stop(msg, " Run brcm_ebm_check(), or set BRCM_STRICT=FALSE for an exploratory run.", call.=FALSE)
        warning(msg, call.=FALSE)}
      else cat(sprintf("EBM reachable in all %d workers\n", length(chk)))}
    on.exit(try(parallel::stopCluster(cl),silent=TRUE),add=TRUE)}
  raw<-list();dg<-list();itr<-list();t0<-Sys.time()
  pr<-.prog_new(length(datasets),sprintf("Benchmark datasets (%d repeats x %d folds)",n_repeats,n_folds))
  ckdir<-file.path(out_dir,"_checkpoints");dir.create(ckdir,recursive=TRUE,showWarnings=FALSE)
  ck_file<-function(x) file.path(ckdir,paste0(x,".rds"))
  mf<-file.path(ckdir,"_manifest.rds")
  cfg<-list(n_repeats=n_repeats,n_folds=n_folds,base_seed=base_seed,
            n_methods=length(registry),datasets=names(datasets))
  if(file.exists(mf)){
    oc<-tryCatch(readRDS(mf),error=function(e)NULL)
    if(!is.null(oc)){
      diffs<-character(0)
      for(k in names(cfg)) if(!identical(oc[[k]],cfg[[k]]))
        diffs<-c(diffs,sprintf("%s: checkpoints used %s, this run uses %s",
                               k,paste(oc[[k]],collapse=","),paste(cfg[[k]],collapse=",")))
      if(length(diffs))
        stop("The existing checkpoints in ",ckdir," do not match this run:\n  - ",
             paste(diffs,collapse="\n  - "),
             "\nUse a different BRCM_OUT, or delete that directory to start fresh.")}
  } else saveRDS(cfg,mf)
  done_before<-length(list.files(ckdir,pattern="\\.rds$"))-1L
  if(done_before>0) cat(sprintf("Resuming: %d of %d datasets already on disk in %s\n",
                                done_before,length(datasets),ckdir))
  for(nm in names(datasets)){
    f<-ck_file(nm)
    if(file.exists(f)){
      r<-tryCatch(readRDS(f),error=function(e)NULL)
      if(!is.null(r)){
        raw[[nm]]<-r$raw;dg[[nm]]<-r$diag
        if(!is.null(r$inter)) itr[[nm]]<-r$inter
        pr<-.prog_tick(pr,sprintf("%s (from checkpoint)",nm));next}}
    r<-run_dataset(datasets[[nm]],registry,n_repeats,n_folds,base_seed,cl,
                   progress=TRUE,progress_label=sprintf("%s repeats",nm))
    if(!is.null(r$cl) && !identical(r$cl,cl)) cl<-r$cl        
    raw[[nm]]<-r$raw;r$diag$Dataset<-nm;dg[[nm]]<-r$diag
    if(!is.null(r$inter)){r$inter$Dataset<-nm;itr[[nm]]<-r$inter}
    tmp<-paste0(f,".tmp")
    saveRDS(list(raw=r$raw,diag=r$diag,inter=r$inter),tmp);file.rename(tmp,f)
    pr<-.prog_tick(pr,sprintf("%s (n=%d, p=%d, prev=%.3f, %d methods)",nm,
      nrow(datasets[[nm]]$X),ncol(datasets[[nm]]$X),mean(datasets[[nm]]$y),length(r$raw)))}
  .prog_done(pr)
  P<-paired_analysis(raw)
  IT<-if(length(itr)) do.call(rbind,itr) else NULL
  isum<-NULL
  if(!is.null(IT)){
    isum<-do.call(rbind,lapply(split(IT,IT$Dataset),function(g){
      structurally_disabled<-all(!isTRUE(any(g$search_enabled)))
      data.frame(Dataset=g$Dataset[1],Outer_fits=nrow(g),
        Any_interaction_pct=round(100*mean(g$n_inter>0L),1),
        Mean_n_inter=round(mean(g$n_inter),3),
        Median_n_inter=stats::median(g$n_inter),
        Mean_candidate_union_pairs=round(mean(g$n_screened_pairs,na.rm=TRUE),1),
        Structurally_disabled=structurally_disabled,stringsAsFactors=FALSE)}))
    rownames(isum)<-NULL}
  res<-list(raw=raw,paired=P,diag=do.call(rbind,dg),inter_fold=IT,inter_summary=isum,
    meta=do.call(rbind,lapply(names(datasets),function(nm){D<-datasets[[nm]];data.frame(Dataset=nm,
      Label=D$label,n=nrow(D$X),p=ncol(D$X),Event_rate=round(mean(D$y),3),
      Missing_pct=round(100*mean(is.na(D$X)),2),
      Predictor_type=if(!is.null(D$ptype))D$ptype else NA,
      Missing_rule=if(!is.null(D$miss))D$miss else NA,
      Outcome_encoding=if(!is.null(D$enc))D$enc else NA,stringsAsFactors=FALSE)})),
    config=list(n_repeats=n_repeats,n_folds=n_folds,oof_K=OOF_K,base_seed=base_seed,
      methods=names(registry),sessionInfo=utils::capture.output(utils::sessionInfo())))
  saveRDS(res,file.path(out_dir,"realdata.rds"))
  dir.create(file.path(out_dir,"tables"),recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(res$meta,file.path(out_dir,"tables","RT0_datasets.csv"),row.names=FALSE)
  write_tables_rd(res,out_dir);write_figures_rd(res,datasets,out_dir)
  cat("Elapsed:",format(round(difftime(Sys.time(),t0),1)),"| output in",normalizePath(out_dir),"\n")
  invisible(res)}

brcm_self_test <- function(verbose=TRUE){
  ok<-function(s,v){if(verbose)cat(sprintf("  %-52s %s\n",s,if(isTRUE(v))"PASS" else "*** FAIL ***"));isTRUE(v)}
  res<-logical(0);ex<-function(z)1/(1+exp(-z))
  set.seed(11);x<-c(1,2,2,2,5,7,7,9)
  res<-c(res,ok("rank map consistent: u_train(x) == u_test(x,x)",max(abs(u_train(x)-u_test(x,x)))<1e-12))
  set.seed(12);n<-900;v1<-rbinom(n,1,.5);v2<-rbinom(n,1,.5);v3<-rnorm(n);v4<-rnorm(n)
  X<-cbind(v1,v2,v3,v4);cm<-c(TRUE,TRUE,FALSE,FALSE)
  y<-rbinom(n,1,ex(-.2+.3*v1+.4*v3+.4*v4+1.6*v3*v4))
  f<-brcm_fit(X,y,cm)
  res<-c(res,ok("interaction selected on a true interaction",length(f$inter)>=1))
  if(length(f$inter)){
    D<-map_main_test(X,f$meta,f$m);it<-f$inter[[1]];A<-.pair_A_new(D,f$block,it$j,it$k)
    Tp<-.pair_T_new(X,f$meta,it$j,it$k,f$m2)-A%*%it$Theta
    res<-c(res,ok("pure interaction orthogonal to [1,main_j,main_k]",max(abs(crossprod(A,Tp)))/nrow(A)<1e-8))
    cp<-brcm_components(f,X)
    eta_c<-cp$intercept+rowSums(cp$main)+rowSums(cp$inter)
    eta_d<-drop(cbind(1,cbind(D,Tp))%*%f$beta)
    res<-c(res,ok("prediction == sum of the plotted components",max(abs(eta_c-eta_d))<1e-10))
    off<-1L+length(f$block);g<-f$beta[off+seq_len(it$ncol)]
    uj<-u_test(X[,it$j],f$meta[[it$j]]$xsort);uk<-u_test(X[,it$k],f$meta[[it$k]]$xsort)
    surf<-drop((.tensor(uj,uk,f$m2)-cbind(1,sweep(bern(uj,f$m),2,f$meta[[it$j]]$ctr),
      sweep(bern(uk,f$m),2,f$meta[[it$k]]$ctr))%*%it$Theta)%*%g)
    res<-c(res,ok("plotted interaction surface == model contribution",max(abs(surf-cp$inter[,1]))<1e-10))
    j<-1;i<-which(f$block==j);bj<-f$beta[1+i];mj<-f$meta[[j]]
    pc<-c(0,bj)-sum(mj$ctr*bj)
    dc<-vapply(c(mj$ref,mj$levels),function(v){z<-as.numeric(mj$levels==v);sum((z-mj$ctr)*bj)},0)
    res<-c(res,ok("plotted categorical bars == model contribution",max(abs(pc-dc))<1e-12))}
  ym<-rbinom(n,1,ex(-.2+.8*v1-.6*v2+1.0*v3));fm<-brcm_fit(X,ym,cm);cm2<-brcm_components(fm,X)
  res<-c(res,ok("continuous components sum-to-zero on training data",max(abs(colMeans(cm2$main[,3:4,drop=FALSE])))<1e-8))
  a<-brcm_fit(X,ym,cm);b<-brcm_fit(X,ym,cm)
  res<-c(res,ok("deterministic (consumes no random numbers)",isTRUE(all.equal(a$beta,b$beta))))
  Xs<-X;Xs[,1]<-1
  res<-c(res,ok("single-level categorical handled",!inherits(try(brcm_fit(Xs,ym,cm),silent=TRUE),"try-error")))
  Xn<-X[1:5,,drop=FALSE];Xn[,1]<-9
  res<-c(res,ok("unseen test category handled",all(is.finite(brcm_predict(fm,Xn)))))
  yn<-rbinom(n,1,0.35);fn<-brcm_fit(X,yn,cm)
  fn<-brcm_fit(X,y,cm,r2_thresh=Inf)
  res<-c(res,ok("null option: r2_thresh=Inf yields zero interactions",fn$n_inter==0L))
  set.seed(31);nw<-700L;pw<-9L
  Xw<-matrix(rnorm(nw*pw),nw,pw);cmw<-rep(FALSE,pw)
  yw<-rbinom(nw,1,ex(0.3*Xw[,1]+0.3*Xw[,2]+1.5*Xw[,1]*Xw[,2]))
  fw<-brcm_fit(Xw,yw,cmw,top_feat=8L)
  res<-c(res,ok("candidate union never exceeds theoretical pair count",
    is.numeric(fw$n_screened_pairs) && fw$n_screened_pairs>0L &&
    fw$n_screened_pairs <= choose(pw,2)))
  fs<-brcm_fit(X[1:250,,drop=FALSE],y[1:250],cm,min_n_inter=300)
  res<-c(res,ok("interaction search disabled below min_n_inter",
    !isTRUE(fs$interaction_search_enabled)))
  set.seed(77);nd<-300L;pd<-12L
  Xd<-cbind(matrix(rbinom(nd*2L,1,.5),nd,2L),matrix(rnorm(nd*(pd-2L)),nd,pd-2L))
  cmd<-c(TRUE,TRUE,rep(FALSE,pd-2L))
  dd<-.dummy(Xd[1:200,,drop=FALSE],Xd[201:nd,,drop=FALSE],cmd)
  cont_ok<-all(vapply(3:pd,function(j)
    any(apply(dd$te,2,function(cc) isTRUE(all.equal(unname(cc),unname(Xd[201:nd,j]))))),TRUE))
  res<-c(res,ok("dummy expansion preserves every continuous predictor",
    cont_ok && !any(duplicated(colnames(dd$tr))) && ncol(dd$tr)==ncol(dd$te)))
  if(verbose)cat(sprintf("  ---- %d/%d passed ----\n",sum(res),length(res)))
  invisible(all(res))}

if(AUTORUN) res <- run_full_realdata()
