.envnum <- function(k,d){v<-Sys.getenv(k);if(nzchar(v))as.numeric(v) else d}
.envchr <- function(k,d){v<-Sys.getenv(k);if(nzchar(v))v else d}
NUM_SIM   <- .envnum("BRCM_NSIM", 1000)                 
AYAR_ns   <- c(100, 500, 1000)                          
AYAR_prev <- c(0.10, 0.25, 0.50)                        
OOF_K     <- .envnum("BRCM_OOFK", 3)                   
BASE_SEED <- .envnum("BRCM_SEED", 1000)
N_CORES   <- .envnum("BRCM_CORES", max(1, parallel::detectCores() - 1))
OUT_DIR   <- .envchr("BRCM_OUT", "BRCMsim")
TEST_FRAC <- 0.30
R2_THRESH <- .envnum("BRCM_R2THRESH", 0.01)             
R2_SENS   <- c(0, 0.005, 0.01, 0.02, 0.05)              
C_SENS    <- c(0.01, 0.025, 0.05, 0.10, 0.20)           
M_SENS    <- c(6, 8, 10, 12)                            
LAMINT_SENS <- c(0.25, 0.5, 1, 2, 4)                    
MININT_SENS <- c(150, 200, 300, 400, 600)               
PROG_CHUNKS <- as.integer(.envnum("BRCM_PROGRESS", 20))
LL_THRESH <- .envnum("BRCM_LLTHRESH", 0.002)            
AUTORUN   <- toupper(.envchr("BRCM_AUTORUN", "TRUE")) == "TRUE"
STRICT_PKGS <- toupper(.envchr("BRCM_STRICT", "TRUE")) == "TRUE"
REQUIRED_PKGS  <- c("mgcv","rpart","randomForest","e1071")
MANUSCRIPT_PKGS <- c("xgboost","lightgbm","catboost")           
OPTIONAL_PKGS  <- c("foreach","doParallel","openxlsx")

.check_pkgs <- function(){
  need <- if (STRICT_PKGS) c(REQUIRED_PKGS, MANUSCRIPT_PKGS) else REQUIRED_PKGS
  miss <- need[!vapply(need, requireNamespace, TRUE, quietly = TRUE)]
  if (length(miss)) {
    msg <- paste0("Missing packages required for this run: ", paste(miss, collapse = ", "))
    if (STRICT_PKGS) stop(msg, ". Install them for the manuscript run, or set BRCM_STRICT=FALSE for a development run that skips missing methods.")
    warning(msg, " -- those methods will be skipped in this development run.")
  }
  invisible(TRUE)
}
has <- function(p) requireNamespace(p, quietly = TRUE)

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
                     sel_crit=c("r2","logloss"),ll_thresh=LL_THRESH,
                     screen_mode=c("nested","outer"),
                     use_rank=TRUE,use_bernstein=TRUE,use_spearman=TRUE,use_inter=TRUE){
  sel_crit<-match.arg(sel_crit);screen_mode<-match.arg(screen_mode)
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
      rv<-c();rh<-c();ll0<-c();ll1<-c();nf<-0L
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
        rv<-c(rv,r_va);rh<-c(rh,drop(Tvp%*%g))
        if(sel_crit=="logloss"){
          e0<-drop(cbind(1,Dv)%*%ft$beta); e1<-e0+drop(Tvp%*%g)
          ll0<-c(ll0,.ll(y[vai],1/(1+exp(-e0)))); ll1<-c(ll1,.ll(y[vai],1/(1+exp(-e1))))}}
      if(!is.null(fold_ok) && nf < (floor(cv_K/2)+1L)) return(-Inf)
      if(sel_crit=="logloss"){
        if(!length(ll0))return(-Inf)
        return(mean(ll0,na.rm=TRUE)-mean(ll1,na.rm=TRUE))}  
      if(length(rv)<10)return(-Inf);1-sum((rv-rh)^2)/sum((rv-mean(rv))^2)},numeric(1))
    thr_sel<-if(sel_crit=="logloss") ll_thresh else r2_thresh   
    keep<-which(scores>thr_sel);keep<-keep[order(scores[keep],decreasing=TRUE)];keep<-head(keep,Kmax)
    sel<-pairs[keep]
    if(length(keep))selinfo<-data.frame(j=vapply(sel,`[`,0L,1),k=vapply(sel,`[`,0L,2),
      oof_score=round(scores[keep],5),score_type=sel_crit,
      n_folds_screened=if(is.null(fold_ok)) cv_K else vapply(keep,function(i)sum(fold_ok[[i]]),0L))}
  Dfull<-D;ipen<-c();inter<-list()
  for(pr in sel){pp<-tryCatch(pair_pure(pr[1],pr[2],X,mm$meta,m,m2,D,mm$block),error=function(e)NULL)
    if(is.null(pp))next
    Dfull<-cbind(Dfull,pp$T);ipen<-c(ipen,rep(lam_int,pp$ncol));inter[[length(inter)+1]]<-pp}
  ff<-.irls(cbind(1,Dfull),y,c(0,lam*penblk,ipen))
  structure(list(beta=ff$beta,lambda=lam,meta=mm$meta,block=mm$block,m=m,m2=m2,use_rank=use_rank,
    inter=inter,selinfo=selinfo,all_oof_score=scores,sel_crit=sel_crit,
    screen_mode=screen_mode,n_screened_pairs=n_screened_pairs,
    interaction_search_enabled=interaction_search_enabled,
    all_oof_R2=if(sel_crit=="r2") scores else NULL,   
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
  attr(p,"fit_diag")<-list(converged=as.integer(f$converged),iter=f$iter,lambda=f$lambda,n_inter=f$n_inter);p}
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

build_registry <- function(){
  reg <- list("B-RCM++"=m_brcm,"LR"=m_lr,"PLR"=m_plr,"DT"=m_dt)
  if(has("mgcv")) reg[["GAM"]] <- m_gam
  if(has("randomForest")) reg[["RF"]] <- m_rf
  if(has("e1071")) reg[["SVM"]] <- m_svm
  if(has("xgboost")) reg[["XGBoost"]] <- m_xgb
  if(has("lightgbm")) reg[["LightGBM"]] <- m_lgbm
  if(has("catboost")) reg[["CatBoost"]] <- m_catb
  reg[["A1-LinearRidge"]]<-m_abl_lin; reg[["A2-RankLinear"]]<-m_abl_rank
  reg[["A3-RankBernstein"]]<-m_abl_bern; reg[["A4-PlusSpearman"]]<-m_abl_spw
  reg }
CALIBRATE <- c("RF","XGBoost","LightGBM","CatBoost")

.expit <- function(z) 1/(1+exp(-z))
.DGP <- list(
  logit_mixed    = list(link="logit", f=function(n){v1<-rbinom(n,1,.5);v2<-rbinom(n,1,.5);v3<-rnorm(n);v4<-rnorm(n)
                          list(X=cbind(v1,v2,v3,v4), lp=0.8*v1-0.6*v2+1.0*v3+0.0*v4)}),
  logit_nl       = list(link="logit", f=function(n){v1<-rbinom(n,1,.5);v2<-rbinom(n,1,.5);v3<-rnorm(n);v4<-rnorm(n)
                          list(X=cbind(v1,v2,v3,v4), lp=0.8*v1-0.6*v2+1.2*tanh(v3))}),
  logit_corr     = list(link="logit", f=function(n){v1<-rbinom(n,1,.5);v2<-rbinom(n,1,.5)
                          z1<-rnorm(n);z2<-0.6*z1+sqrt(1-0.36)*rnorm(n)
                          list(X=cbind(v1,v2,z1,z2), lp=0.7*v1-0.5*v2+0.9*z1+0.5*z2)}),
  lpm_mixed      = list(link="lpm",   f=function(n){v1<-rbinom(n,1,.5);v2<-rbinom(n,1,.5);v3<-rnorm(n);v4<-rnorm(n)
                          list(X=cbind(v1,v2,v3,v4), lp=0.10*v1-0.08*v2+0.12*v3+0.0*v4)}),
  logit_interact = list(link="logit", f=function(n){v1<-rbinom(n,1,.5);v2<-rbinom(n,1,.5);v3<-rnorm(n);v4<-rnorm(n)
                          list(X=cbind(v1,v2,v3,v4), lp=0.3*v1+0.4*v3+0.4*v4+1.6*v3*v4)}),
  logit_pure_int = list(link="logit", f=function(n){v1<-rbinom(n,1,.5);v2<-rbinom(n,1,.5);v3<-rnorm(n);v4<-rnorm(n)
                          list(X=cbind(v1,v2,v3,v4), lp=1.6*v3*v4)}), 
  null_mixed     = list(link="logit", f=function(n){v1<-rbinom(n,1,.5);v2<-rbinom(n,1,.5);v3<-rnorm(n);v4<-rnorm(n)
                          list(X=cbind(v1,v2,v3,v4), lp=rep(0,n))}),
  logit_ushape   = list(link="logit", f=function(n){v1<-rbinom(n,1,.5);v2<-rbinom(n,1,.5);v3<-rnorm(n);v4<-rnorm(n)
                          list(X=cbind(v1,v2,v3,v4), lp=0.5*v1-0.4*v2+1.4*(v3^2-1)+0.3*v4)}))

.augment_dgp <- function(base, p){
  stopifnot(p > 4L)
  bf <- .DGP[[base]]$f; lk <- .DGP[[base]]$link
  n_noise <- as.integer(p - 4L)
  nb <- as.integer(ceiling(n_noise/2)); nc <- as.integer(n_noise - nb)
  list(link = lk, p = as.integer(p), n_noise_bin = nb, n_noise_cont = nc,
       f = function(n){
         d <- bf(n)                                  
         Z <- cbind(if(nb) matrix(rbinom(n*nb, 1, .5), n, nb) else NULL,
                    if(nc) matrix(rnorm(n*nc), n, nc) else NULL)
         list(X = cbind(d$X, Z), lp = d$lp)})}        

.DGP_P <- c(10L, 20L)                                  
for(.b in names(.DGP)) for(.pp in .DGP_P)
  .DGP[[sprintf("%s_p%d", .b, .pp)]] <- .augment_dgp(.b, .pp)
rm(.b, .pp)

CAT_MASK <- c(TRUE,TRUE,FALSE,FALSE)              
.DGP_CAT <- list()                                
dgp_cat_mask <- function(dgp){
  cm <- .DGP_CAT[[dgp]]
  if(!is.null(cm)) return(cm)
  d <- .DGP[[dgp]]
  if(!is.null(d$n_noise_bin))
    return(c(CAT_MASK, rep(TRUE, d$n_noise_bin), rep(FALSE, d$n_noise_cont)))
  if(!is.null(d$p)) return(rep(FALSE, d$p))          
  CAT_MASK}

.register_highdim <- function(p){
  nm <- paste0("hd",p)
  .DGP[[paste0(nm,"_null")]] <<- list(link="logit", p=p,
    f=function(n){X<-matrix(rnorm(n*p),n,p); list(X=X, lp=rep(0,n))})
  .DGP[[paste0(nm,"_interact")]] <<- list(link="logit", p=p,
    f=function(n){X<-matrix(rnorm(n*p),n,p)
      list(X=X, lp=0.4*X[,1]+0.4*X[,2]+1.6*X[,1]*X[,2])})
  .DGP[[paste0(nm,"_pure_int")]] <<- list(link="logit", p=p,
    f=function(n){X<-matrix(rnorm(n*p),n,p); list(X=X, lp=1.6*X[,1]*X[,2])})
  invisible(NULL)}
for(.p in c(10,20,50)) .register_highdim(.p)
rm(.p)

DGP_LABEL <- c(logit_mixed="Logit (mixed)",logit_nl="Logit (nonlinear)",logit_corr="Logit (correlated)",
  lpm_mixed="LPM (mixed)",logit_interact="Logit (interaction)",logit_pure_int="Logit (pure interaction)",
  null_mixed="Null",logit_ushape="Logit (U-shaped)")
DGP_LABEL <- c(DGP_LABEL, unlist(lapply(.DGP_P, function(pp)
  setNames(sprintf("%s, p=%d", DGP_LABEL, pp),
           sprintf("%s_p%d", names(DGP_LABEL), pp)))))
DGP_P_OF <- vapply(names(DGP_LABEL), function(d){
  pp <- .DGP[[d]]$p; if(is.null(pp)) 4L else as.integer(pp)}, 0L)
HD_LABEL <- unlist(lapply(c(10,20,50), function(p) setNames(
  c(sprintf("HD p=%d (null)",p),sprintf("HD p=%d (interaction)",p),sprintf("HD p=%d (pure interaction)",p)),
  paste0("hd",p,c("_null","_interact","_pure_int")))))
BASE_DGPS <- names(DGP_LABEL)          

.INTERCEPT_CACHE <- new.env(parent=emptyenv())
solve_intercept <- function(dgp,target,pool=200000){
  key<-paste0(dgp,"_",target); if(!is.null(.INTERCEPT_CACHE[[key]])) return(.INTERCEPT_CACHE[[key]])
  old <- if(exists(".Random.seed",envir=.GlobalEnv)) get(".Random.seed",envir=.GlobalEnv) else NULL
  set.seed(12345)
  reg<-.DGP[[dgp]];lp<-reg$f(pool)$lp
  g<-function(b0) if(reg$link=="logit") mean(.expit(b0+lp)) else mean(pmin(pmax(b0+lp,0.001),0.999))
  lo<--20;hi<-20
  for(i in 1:200){mid<-(lo+hi)/2;if(g(mid)<target)lo<-mid else hi<-mid}
  b0<-(lo+hi)/2
  if(!is.null(old)) assign(".Random.seed",old,envir=.GlobalEnv)
  .INTERCEPT_CACHE[[key]]<-b0;b0}

generate_dataset <- function(dgp,n,target){
  reg<-.DGP[[dgp]];d<-reg$f(n);b0<-solve_intercept(dgp,target)
  pr<-if(reg$link=="logit") .expit(b0+d$lp) else pmin(pmax(b0+d$lp,0.001),0.999)
  list(X=d$X,y=rbinom(n,1,pr))}

split_stratified <- function(y,test_frac=TEST_FRAC){
  te<-integer(0)
  for(cl in unique(y)){i<-which(y==cl);k<-max(1,round(length(i)*test_frac));te<-c(te,sample(i,k))}
  list(train=setdiff(seq_along(y),te),test=sort(te))}

oof_probs <- function(fitter,Xtr,ytr,cm,K=OOF_K){
  fo<-.strat_folds(ytr,K);p<-numeric(length(ytr))
  for(f in 1:K){tri<-which(fo!=f);vai<-which(fo==f)
    if(length(unique(ytr[tri]))<2){p[vai]<-mean(ytr[tri]);next}
    pv<-tryCatch(fitter(Xtr[tri,,drop=FALSE],ytr[tri],Xtr[vai,,drop=FALSE],cm),error=function(e)e)
    if(inherits(pv,"error")||length(pv)!=length(vai)||any(!is.finite(pv)))
      stop("OOF inner fit failed")
    p[vai]<-pv}
  p}

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

run_one_rep <- function(dgp,n,target,rep_id,registry,base_seed=BASE_SEED){
  set.seed(base_seed+rep_id)
  d<-generate_dataset(dgp,n,target)
  cmask<-dgp_cat_mask(dgp)              
  sp<-split_stratified(d$y)
  Xtr<-d$X[sp$train,,drop=FALSE];ytr<-d$y[sp$train]
  Xte<-d$X[sp$test,,drop=FALSE];yte<-d$y[sp$test]
  bad<-length(unique(ytr))<2||length(unique(yte))<2
  out<-list();diag<-list();extra<-list()
  midx<-0L
  for(nm in names(registry)){
    midx<-midx+1L
    if(bad){diag[[nm]]<-c(ok=0,reason=1,secs=NA,converged=NA,iter=NA,lambda=NA,n_inter=NA);next}
    set.seed(base_seed*7L + rep_id*1000L + .name_seed(nm))
    fitter<-registry[[nm]]
    t0<-Sys.time()
    res<-tryCatch({
      p_oof<-oof_probs(fitter,Xtr,ytr,cmask)
      thr<-.youden(ytr,p_oof)
      p_te<-fitter(Xtr,ytr,Xte,cmask)
      if(any(!is.finite(p_te))) stop("non-finite predictions")
      list(main=evaluate(yte,p_te,thr),p_oof=p_oof,p_te=p_te,diag=attr(p_te,"fit_diag"))
      },error=function(e){attr(e,"msg")<-conditionMessage(e);e})
    secs<-as.numeric(difftime(Sys.time(),t0,units="secs"))
    if(inherits(res,"error")){diag[[nm]]<-c(ok=0,reason=2,secs=secs,converged=NA,iter=NA,lambda=NA,n_inter=NA);next}
    out[[nm]]<-res$main
    fd<-res$diag
    diag[[nm]]<-c(ok=1,reason=0,secs=secs,
      converged=if(is.null(fd))NA else fd$converged,iter=if(is.null(fd))NA else fd$iter,
      lambda=if(is.null(fd))NA else fd$lambda,n_inter=if(is.null(fd))NA else fd$n_inter)
    if(nm %in% CALIBRATE){
      cnm<-paste0(nm,"+Platt")
      cal<-tryCatch({
        cf<-.platt_fit(ytr,res$p_oof)                        
        thr_c<-.youden(ytr,.platt_oof(ytr,res$p_oof))       
        evaluate(yte,.platt_predict(cf,res$p_te),thr_c)},error=function(e)e)
      if(inherits(cal,"error")){diag[[cnm]]<-c(ok=0,reason=3,secs=NA,converged=NA,iter=NA,lambda=NA,n_inter=NA)}
      else {out[[cnm]]<-cal;diag[[cnm]]<-c(ok=1,reason=0,secs=NA,converged=NA,iter=NA,lambda=NA,n_inter=NA)}}
  }
  list(metrics=out,diag=diag,bad_split=bad,rep_id=rep_id)
}

run_scenario <- function(dgp,n,target,num_sim,registry,base_seed=BASE_SEED,cl=NULL,
                        progress=TRUE,progress_chunks=PROG_CHUNKS,progress_label="",
                        n_cores=1L){
  reps<-seq_len(num_sim)
  runner<-function(r) run_one_rep(dgp,n,target,r,registry,base_seed)
  {
    quiet<-!(isTRUE(progress) && progress_chunks >= 1L)
    nchunk<-if(quiet) 1L else max(1L, min(as.integer(progress_chunks), num_sim))
    bounds<-unique(as.integer(round(seq(0, num_sim, length.out = nchunk + 1L))))
    lbl<-if(nzchar(progress_label)) progress_label else sprintf("%s reps",dgp)
    t0<-Sys.time(); res<-vector("list", num_sim)
    .newcl<-function(){
      k<-tryCatch(parallel::makeCluster(n_cores),error=function(e)NULL)
      if(is.null(k)) return(NULL)
      ok<-tryCatch({parallel::clusterExport(k,varlist=ls(envir=.GlobalEnv,all.names=TRUE),
                                            envir=.GlobalEnv)
                    parallel::clusterEvalQ(k,{suppressMessages({library(stats);library(utils)})})
                    TRUE},error=function(e)FALSE)
      if(!ok){try(parallel::stopCluster(k),silent=TRUE);return(NULL)}
      k}
    for(b in seq_len(length(bounds)-1L)){
      idx<-(bounds[b]+1L):bounds[b+1L]
      part<-tryCatch(if(!is.null(cl)) parallel::parLapply(cl,idx,runner) else lapply(idx,runner),
                     error=function(e)NULL)
      if(is.null(part)){
        cat(sprintf("\n    [!] chunk %d of %s lost a worker; rebuilding and isolating\n",
                    b, if(nzchar(progress_label)) progress_label else dgp));flush.console()
        if(!is.null(cl)) try(parallel::stopCluster(cl),silent=TRUE)
        cl<-.newcl()
        part<-vector("list",length(idx))
        for(u in seq_along(idx)){
          one<-tryCatch({
            k1<-parallel::makeCluster(1L)
            on.exit(try(parallel::stopCluster(k1),silent=TRUE),add=TRUE)
            parallel::clusterExport(k1,varlist=ls(envir=.GlobalEnv,all.names=TRUE),envir=.GlobalEnv)
            v<-parallel::parLapply(k1,idx[u],runner)[[1]]
            try(parallel::stopCluster(k1),silent=TRUE);v},
            error=function(e)NULL)
          part[[u]]<-one}
        nlost<-sum(vapply(part,is.null,TRUE))
        if(nlost) cat(sprintf("    [!] %d replication(s) in this chunk could not be completed\n",nlost))}
      res[idx]<-part
      dn<-bounds[b+1L]
      el<-as.numeric(difftime(Sys.time(),t0,units="secs"))
      eta<-if(dn>0L) el/dn*(num_sim-dn) else NA_real_
      if(!quiet){
        cat(sprintf("\r    %s %d/%d (%5.1f%%) | elapsed %s | eta %s   ",
                    lbl, dn, num_sim, 100*dn/num_sim, .fmt_hms(el), .fmt_hms(eta)))
        flush.console()}}
    if(!quiet) cat("\n")
  }
  keep<-!vapply(res,is.null,TRUE)
  if(!all(keep)) res<-res[keep]
  allm<-unique(unlist(lapply(res,function(z)names(z$metrics))))
  raw<-list()
  for(nm in allm){
    ids<-integer(0);rows<-list()
    for(z in res) if(!is.null(z$metrics[[nm]])){ids<-c(ids,z$rep_id);rows[[length(rows)+1]]<-z$metrics[[nm]]}
    if(!length(rows))next
    M<-do.call(rbind,rows);rownames(M)<-as.character(ids)     
    raw[[nm]]<-M}
  getd<-function(z,nm,f){v<-z$diag[[nm]];if(is.null(v)||is.na(v[f]))NA_real_ else as.numeric(v[f])}
  diag_methods<-unique(c(names(registry),unlist(lapply(res,function(z)names(z$diag)))))
  dg<-do.call(rbind,lapply(diag_methods,function(nm){
    att<-length(res);okn<-sum(vapply(res,function(z)isTRUE(!is.null(z$metrics[[nm]])),TRUE))
    badsp<-sum(vapply(res,function(z)isTRUE(z$bad_split),TRUE))
    secs<-vapply(res,getd,0,nm=nm,f="secs")
    data.frame(Method=nm,Attempted=att,Successful=okn,Bad_split=badsp,
      Failure_rate=round(100*(att-okn)/att,2),
      Median_secs=round(.safe_median(secs),4),
      Converged_pct=round(100*.safe_mean(vapply(res,getd,0,nm=nm,f="converged")),1),
      Median_iter=round(.safe_median(vapply(res,getd,0,nm=nm,f="iter")),1),
      Median_lambda=round(.safe_median(vapply(res,getd,0,nm=nm,f="lambda")),3),
      Mean_n_inter=round(.safe_mean(vapply(res,getd,0,nm=nm,f="n_inter")),2),
      stringsAsFactors=FALSE)}))
  list(raw=raw,diag=dg,cl=cl)
}

.rank_biserial <- function(d){d<-d[is.finite(d)];d<-d[d!=0];if(length(d)<2)return(NA_real_)
  r<-rank(abs(d));(sum(r[d>0])-sum(r[d<0]))/sum(r)}
LOWER_BETTER <- c(Test_AUC=FALSE,Test_PRAUC=FALSE,Test_Brier=TRUE,Test_LogLoss=TRUE,
  Test_ECE=TRUE,Test_ECEadapt=TRUE,Test_ICI=TRUE)

paired_analysis <- function(raw_all,reference="B-RCM++",metrics=names(LOWER_BETTER)){
  rows<-list();k<-0
  for(scn in names(raw_all)){
    R<-raw_all[[scn]];ref<-R[[reference]];if(is.null(ref))next
    for(cmp in setdiff(names(R),reference)){
      oth<-R[[cmp]];if(is.null(oth))next
      ids<-intersect(rownames(ref),rownames(oth))         
      if(length(ids)<3)next
      for(met in metrics){
        if(!(met %in% colnames(ref))||!(met %in% colnames(oth)))next
        d<-ref[ids,met]-oth[ids,met];d<-d[is.finite(d)]
        if(length(d)<3)next
        dd<-if(LOWER_BETTER[[met]]) -d else d               
        wt<-tryCatch(stats::wilcox.test(dd,conf.int=TRUE,exact=FALSE),error=function(e)NULL)
        k<-k+1
        rows[[k]]<-data.frame(Scenario=scn,Comparator=cmp,Metric=met,
          n_common=length(ids),n_used=length(d),
          HL=if(!is.null(wt))unname(wt$estimate) else stats::median(dd),
          CI_lo=if(!is.null(wt))wt$conf.int[1] else NA_real_,
          CI_hi=if(!is.null(wt))wt$conf.int[2] else NA_real_,
          rank_biserial=.rank_biserial(dd),
          p_raw=if(!is.null(wt))wt$p.value else NA_real_,stringsAsFactors=FALSE)}}}
  D<-do.call(rbind,rows);if(is.null(D))return(NULL)
  D$p_BH<-NA_real_
  for(cmp in unique(D$Comparator)) for(met in unique(D$Metric)){
    i<-which(D$Comparator==cmp&D$Metric==met);if(length(i))D$p_BH[i]<-stats::p.adjust(D$p_raw[i],method="BH")}
  D$p_BH_global  <- stats::p.adjust(D$p_raw, method="BH")
  D$p_holm_global<- stats::p.adjust(D$p_raw, method="holm")
  D}

scenario_summary <- function(P){
  do.call(rbind,lapply(split(P,list(P$Comparator,P$Metric),drop=TRUE),function(g)
    data.frame(Comparator=g$Comparator[1],Metric=g$Metric[1],n_scenarios=nrow(g),
      median_HL=stats::median(g$HL,na.rm=TRUE),
      IQR_HL=paste0(round(stats::quantile(g$HL,.25,na.rm=TRUE),4),", ",round(stats::quantile(g$HL,.75,na.rm=TRUE),4)),
      pct_favour_ref=round(100*mean(g$HL>0,na.rm=TRUE),1),
      pct_sig_favour=round(100*mean(g$HL>0&g$p_BH<0.05,na.rm=TRUE),1),
      pct_sig_against=round(100*mean(g$HL<0&g$p_BH<0.05,na.rm=TRUE),1),
      pct_sig_favour_global=round(100*mean(g$HL>0&g$p_BH_global<0.05,na.rm=TRUE),1),
      pct_sig_against_global=round(100*mean(g$HL<0&g$p_BH_global<0.05,na.rm=TRUE),1),
      pct_sig_favour_holm=round(100*mean(g$HL>0&g$p_holm_global<0.05,na.rm=TRUE),1),
      stringsAsFactors=FALSE)))}

.agg <- function(raw_all){
  rows<-list();k<-0
  for(scn in names(raw_all)) for(nm in names(raw_all[[scn]])){
    M<-raw_all[[scn]][[nm]]
    for(met in colnames(M)){v<-M[,met];v<-v[is.finite(v)];if(!length(v))next
      k<-k+1;rows[[k]]<-data.frame(Scenario=scn,Model=nm,Metric=met,N=length(v),
        Mean=mean(v),MC_SE=stats::sd(v)/sqrt(length(v)),SD=stats::sd(v),
        Median=stats::median(v),Q1=stats::quantile(v,.25,names=FALSE),
        Q3=stats::quantile(v,.75,names=FALSE),stringsAsFactors=FALSE)}}
  do.call(rbind,rows)}
.parse_scn <- function(s){m<-regmatches(s,regexec("^(.*)_n([0-9]+)_p([0-9]+)$",s))[[1]]
  list(dgp=m[2],n=as.integer(m[3]),prev=as.integer(m[4]))}

write_tables <- function(sim,out_dir=OUT_DIR){
  td<-file.path(out_dir,"tables");dir.create(td,recursive=TRUE,showWarnings=FALSE)
  S<-.agg(sim$raw);pp<-t(vapply(S$Scenario,function(s){z<-.parse_scn(s);c(z$dgp,z$n,z$prev)},character(3)))
  S$DGP<-DGP_LABEL[pp[,1]];S$n<-as.integer(pp[,2]);S$prev<-as.integer(pp[,3])
  utils::write.csv(S,file.path(td,"DATA_summary.csv"),row.names=FALSE)
  models<-unique(S$Model)
  o<-do.call(rbind,lapply(models,function(m){r<-S[S$Model==m,]
    g<-function(met)mean(r$Mean[r$Metric==met],na.rm=TRUE)
    gmed<-function(met)stats::median(r$Median[r$Metric==met],na.rm=TRUE)
    data.frame(Model=m,AUC=g("Test_AUC"),PR_AUC=g("Test_PRAUC"),Brier=g("Test_Brier"),
      LogLoss=g("Test_LogLoss"),ICI=g("Test_ICI"),ECEadapt=g("Test_ECEadapt"),
      ECE=g("Test_ECE"),CalSlope=gmed("Test_CalSlope"),
      BalAcc_Y=g("Test_BalAcc_Y"),stringsAsFactors=FALSE)}))
  au<-S[S$Metric=="Test_AUC",];w<-stats::reshape(au[,c("Scenario","Model","Mean")],idvar="Scenario",
    timevar="Model",direction="wide");rk<-t(apply(-as.matrix(w[,-1]),1,rank,na.last="keep"))
  colnames(rk)<-sub("^Mean\\.","",colnames(w)[-1])
  o$MeanRank_AUC<-round(colMeans(rk,na.rm=TRUE)[o$Model],2)
  utils::write.csv(o,file.path(td,"T1_overall.csv"),row.names=FALSE)
  fmt<-function(md,q1,q3,pct)if(is.na(md))"--" else if(pct)sprintf("%.1f (%.1f)",100*md,100*(q3-q1)) else sprintf("%.3f (%.3f)",md,q3-q1)
  METS<-list(c("Test_AUC","AUC","1"),c("Test_PRAUC","PR-AUC","1"),c("Test_Brier","Brier","0"),
    c("Test_LogLoss","Log-loss","0"),c("Test_ICI","ICI","0"),
    c("Test_ECEadapt","ECE (adaptive)","0"),c("Test_ECE","ECE (equal-width)","0"),
    c("Test_CalSlope","Cal. slope","0"),
    c("Test_CalIntercept","Cal. intercept","0"),c("Test_BalAcc_Y","Bal. acc.","1"),
    c("Test_Sens_Y","Sens.","1"),c("Test_Spec_Y","Spec.","1"),c("Test_F1_Y","F1","1"),c("Test_MCC_Y","MCC","0"))
  ti<-2
  for(dg in BASE_DGPS) for(pv in c(10,25,50)){
    blk<-list()
    for(nn in AYAR_ns){scn<-paste0(dg,"_n",nn,"_p",pv);sub<-S[S$Scenario==scn,];if(!nrow(sub))next
      for(me in METS){r<-data.frame(n=nn,Metric=me[2],stringsAsFactors=FALSE)
        for(m in models){x<-sub[sub$Model==m&sub$Metric==me[1],]
          r[[m]]<-if(nrow(x))fmt(x$Median[1],x$Q1[1],x$Q3[1],me[3]=="1") else "--"}
        blk[[length(blk)+1]]<-r}}
    if(length(blk)){utils::write.csv(do.call(rbind,blk),
      file.path(td,sprintf("T%02d_detail_%s_p%02d.csv",ti,dg,pv)),row.names=FALSE);ti<-ti+1}}
  if(!is.null(sim$paired)){
    P<-sim$paired
    P$cell<-sprintf("%+.4f%s (%+.4f, %+.4f) [%+.2f]",P$HL,ifelse(!is.na(P$p_BH_global)&P$p_BH_global<0.05,"*",""),
      P$CI_lo,P$CI_hi,P$rank_biserial)
    utils::write.csv(P,file.path(td,"DATA_paired_scenario.csv"),row.names=FALSE)
    utils::write.csv(scenario_summary(P),file.path(td,"T_paired_summary.csv"),row.names=FALSE)}
  utils::write.csv(sim$diag,file.path(td,"T_diagnostics.csv"),row.names=FALSE)
  if(!is.null(sim$diag)){D<-sim$diag
    rt<-do.call(rbind,lapply(unique(D$Method),function(m){r<-D[D$Method==m,]
      data.frame(Method=m,
        Median_secs_per_eval=round(.safe_median(r$Median_secs),4),
        Mean_failure_rate=round(.safe_mean(r$Failure_rate),2),
        Converged_pct=round(.safe_mean(r$Converged_pct),1),
        Median_lambda=round(.safe_median(r$Median_lambda),3),
        Mean_n_inter=round(.safe_mean(r$Mean_n_inter),2),stringsAsFactors=FALSE)}))
    utils::write.csv(rt,file.path(td,"T_runtime.csv"),row.names=FALSE)}
  if(!is.null(sim$inter_diag)) utils::write.csv(sim$inter_diag,file.path(td,"T_interaction_diagnostics.csv"),row.names=FALSE)
  if(!is.null(sim$highdim)) utils::write.csv(sim$highdim,file.path(td,"T_highdim_screening.csv"),row.names=FALSE)
  if(has("openxlsx")){
    wb<-openxlsx::createWorkbook();add<-function(n,d){openxlsx::addWorksheet(wb,n);openxlsx::writeData(wb,n,d)}
    add("Overall",o);add("Summary_long",S)
    if(!is.null(sim$paired)){add("Paired_scenario",sim$paired);add("Paired_summary",scenario_summary(sim$paired))}
    add("Diagnostics",sim$diag)
    openxlsx::saveWorkbook(wb,file.path(out_dir,"Paper_Tables.xlsx"),overwrite=TRUE)}
  invisible(td)}

draw_screening_figure <- function(H, file){
  H$fam <- gsub(".*\\((.*)\\).*", "\\1", H$Process)
  grDevices::png(file, width=3000, height=1500, res=300)
  op <- graphics::par(mfrow=c(1,2), mar=c(4.2,4.4,2.6,1), mgp=c(2.5,0.8,0))
  on.exit({graphics::par(op); grDevices::dev.off()}, add=TRUE)
  nl <- H[H$fam=="null",]; ns <- sort(unique(nl$n)); cols <- c("#1f77b4","#d62728")
  plot(range(H$p), c(0,100), type="n", xlab="Number of predictors (p)",
       ylab="Fits with at least one false surface (%)", main="Null process: false selection")
  for(i in seq_along(ns)){ x <- nl[nl$n==ns[i],]; x <- x[order(x$p),]
    lines(x$p, x$FWER_false_selection_pct, type="b", pch=16, lwd=2, col=cols[i]) }
  legend("topright", legend=paste0("n = ", ns), col=cols, lwd=2, pch=16, bty="n", cex=0.85)
  grid(col="grey85")
  it <- H[H$fam!="null" & !is.na(H$Correct_pair_pct),]
  fams <- unique(it$fam); cols2 <- c("#2ca02c","#9467bd"); lt <- c(1,2)
  plot(range(H$p), c(0,100), type="n", xlab="Number of predictors (p)",
       ylab="Replications recovering the generating pair (%)", main="Interaction processes: recovery")
  leg <- c(); lc <- c(); ll <- c()
  for(i in seq_along(fams)) for(k in seq_along(ns)){
    x <- it[it$fam==fams[i] & it$n==ns[k],]; x <- x[order(x$p),]
    if(!nrow(x)) next
    lines(x$p, x$Correct_pair_pct, type="b", pch=17, lwd=2, col=cols2[i], lty=lt[k])
    leg <- c(leg, sprintf("%s, n = %d", fams[i], ns[k])); lc <- c(lc, cols2[i]); ll <- c(ll, lt[k]) }
  legend("topright", legend=leg, col=lc, lty=ll, lwd=2, pch=17, bty="n", cex=0.72)
  grid(col="grey85")
  invisible(file)}

write_figures <- function(sim,out_dir=OUT_DIR){
  fd<-file.path(out_dir,"figures");dir.create(fd,recursive=TRUE,showWarnings=FALSE)
  S<-.agg(sim$raw);pp<-t(vapply(S$Scenario,function(s){z<-.parse_scn(s);c(z$dgp,z$n,z$prev)},character(3)))
  S$DGP<-DGP_LABEL[pp[,1]];S$n<-as.integer(pp[,2])
  core<-intersect(c("B-RCM++","LR","PLR","GAM","DT","RF","SVM","XGBoost","LightGBM","CatBoost"),unique(S$Model))
  cols<-grDevices::hcl.colors(length(core),"Dark 3")
  png(file.path(fd,"F1_AUC_by_n.png"),width=2800,height=1800,res=300)
  a<-S[S$Metric=="Test_AUC",];ylim<-range(a$Mean,na.rm=TRUE)
  plot(NA,xlim=range(AYAR_ns),ylim=ylim,log="x",xlab="Sample size (n)",ylab="Mean test AUC",
    main="Discrimination by sample size",xaxt="n");axis(1,at=AYAR_ns,labels=AYAR_ns)
  for(i in seq_along(core)){v<-vapply(AYAR_ns,function(nn)mean(a$Mean[a$Model==core[i]&a$n==nn],na.rm=TRUE),0)
    lines(AYAR_ns,v,col=cols[i],lwd=if(core[i]=="B-RCM++")3 else 1.5,type="b",pch=16)}
  legend("bottomright",legend=core,col=cols,lwd=2,cex=.7,bty="n");dev.off()
  png(file.path(fd,"F2_LogLoss_by_DGP.png"),width=3000,height=2100,res=300)
  l<-S[S$Metric=="Test_LogLoss",];dl<-unique(l$DGP)
  M<-vapply(core,function(m)vapply(dl,function(d)mean(l$Mean[l$Model==m&l$DGP==d],na.rm=TRUE),0),numeric(length(dl)))
  op<-par(mar=c(13,4.4,3,1),mgp=c(2.8,0.8,0))
  barplot(t(M),beside=TRUE,col=cols,las=2,cex.names=0.72,ylab="Mean test log-loss",
    main="Probabilistic performance by process (lower is better)")
  legend("topright",legend=core,fill=cols,cex=.6,bty="n");par(op);dev.off()
  if(!is.null(sim$paired)){
    P<-sim$paired[sim$paired$Metric=="Test_AUC",]
    if(nrow(P)){png(file.path(fd,"F3_HL_forest_AUC.png"),width=2400,height=3000,res=300)
      P<-P[order(P$HL),];op<-par(mar=c(5,11,3,2))
      plot(P$HL,seq_len(nrow(P)),pch=16,cex=.5,yaxt="n",xaxt="n",ylab="",
        xlab="Hodges\u2013Lehmann difference in AUC (positive favours B-RCM++)",
        main="Scenario-level paired effects",xlim=range(c(P$CI_lo,P$CI_hi),na.rm=TRUE))
      axis_minus(1)
      segments(P$CI_lo,seq_len(nrow(P)),P$CI_hi,seq_len(nrow(P)),col="grey40")
      abline(v=0,lty=2,col="red");par(op);dev.off()}}

  if(!is.null(sim$highdim) && nrow(sim$highdim))
    draw_screening_figure(sim$highdim, file.path(fd, "F6_screening_by_dimension.png"))
  set.seed(BASE_SEED)
  d<-generate_dataset("logit_nl",1000,0.50);f<-brcm_fit(d$X,d$y,CAT_MASK)
  brcm_plot_shapes(f,c("V1 (binary)","V2 (binary)","V3 (continuous)","V4 (null)"),
    file=file.path(fd,"F4_component_functions.png"))
  set.seed(BASE_SEED)
  d2<-generate_dataset("logit_interact",1000,0.50);f2<-brcm_fit(d2$X,d2$y,CAT_MASK)
  if(length(f2$inter)) brcm_plot_interaction(f2,1,paste0("V",1:4),file=file.path(fd,"F5_interaction_surface.png"))
  invisible(fd)}

highdim_analysis <- function(num_sim=200, ps=c(10,20,50), ns=c(500,1000), target=0.30,
                             thresholds=c(0.01), top_feats=8L,
                             base_seed=BASE_SEED, verbose=TRUE, n=NULL){
  if(!is.null(n)) ns <- n                      
  is_pair12 <- function(f){
    if(f$n_inter < 1) return(FALSE)
    any(vapply(f$inter, function(it) all(sort(c(it$j,it$k)) == c(1,2)), TRUE))}
  out <- list(); k <- 0L
  .n_cells <- length(ns) * sum(vapply(ps, function(pp){
    tfs0 <- if(is.null(top_feats)) c(8L, min(pp, 16L)) else top_feats
    length(unique(pmin(tfs0, pp))) * length(thresholds) * 3L}, 0L))
  pr_hd <- .prog_new(.n_cells, sprintf("High-dimensional screening (%d reps/cell)", num_sim),
                     quiet = !verbose)
  for(nn in ns) for(p in ps){
    tfs <- if(is.null(top_feats)) c(8L, min(p, 16L)) else top_feats
    tfs <- unique(pmin(tfs, p))
    for(tf in tfs) for(th in thresholds){
      for(fam in c("null","interact","pure_int")){
        dg <- sprintf("hd%d_%s", p, fam)
        cm <- dgp_cat_mask(dg)
        nsucc <- 0L; anysel <- 0L; anyfalse <- 0L; correct <- 0L
        ninter <- c(); nfalse_v <- c(); auc <- c(); ll <- c()
        nunion <- c()                                  
        pr_rep <- .prog_new(num_sim, sprintf("  p=%d %s", p, fam), quiet = TRUE)
        .rep_step <- max(1L, as.integer(round(num_sim/4)))   
        for(r in 1:num_sim){
          if(verbose && r %% .rep_step == 0L)
            cat(sprintf("    ... p=%d n=%d %s: rep %d/%d (%.0f%%)\n", p, nn, fam, r, num_sim, 100*r/num_sim))
          set.seed(base_seed + 5000L + r)
          d  <- generate_dataset(dg, nn, target)
          sp <- split_stratified(d$y)
          f <- tryCatch(brcm_fit(d$X[sp$train,,drop=FALSE], d$y[sp$train], cm,
                                 r2_thresh=th, top_feat=tf), error=function(e) NULL)
          if(is.null(f) || !isTRUE(f$converged)) next
          nsucc  <- nsucc + 1L
          has_true <- if(fam == "null") 0L else as.integer(is_pair12(f))
          nfalse   <- f$n_inter - has_true
          anysel   <- anysel + (f$n_inter > 0)          
          anyfalse <- anyfalse + (nfalse > 0)           
          nfalse_v <- c(nfalse_v, nfalse)
          correct<- correct + is_pair12(f)
          ninter <- c(ninter, f$n_inter)
          nunion <- c(nunion, f$n_screened_pairs)
          pp <- tryCatch(brcm_predict(f, d$X[sp$test,,drop=FALSE]), error=function(e) NULL)
          if(!is.null(pp)){ yt <- d$y[sp$test]; auc <- c(auc, .auc(yt,pp)); ll <- c(ll, .ll(yt,pp)) }
        }
        k <- k + 1L
        out[[k]] <- data.frame(
          p = p, n = nn, p_over_n = round(p/nn, 4),
          n_possible_pairs = p*(p-1L)/2L,          
          screen_size = tf,
          pairs_per_fold = choose(tf, 2L),         
          Mean_candidate_union_pairs = round(.safe_mean(nunion), 1),  
          r2_thresh = th,
          Process = HD_LABEL[[dg]], Attempted = num_sim, Successful = nsucc,
          Pct_any_selection      = if(nsucc) round(100*anysel/nsucc, 1) else NA_real_,
          FWER_false_selection_pct = if(nsucc) round(100*anyfalse/nsucc, 1) else NA_real_,
          Correct_pair_pct       = if(nsucc && fam != "null") round(100*correct/nsucc, 1) else NA_real_,
          Mean_n_surfaces        = if(nsucc) round(.safe_mean(ninter), 3) else NA_real_,
          True_pairs_available   = if(fam == "null") 0L else 1L,
          Mean_false_surfaces    = if(nsucc) round(.safe_mean(nfalse_v), 3) else NA_real_,
          Mean_AUC               = round(.safe_mean(auc), 4),
          Mean_LogLoss           = round(.safe_mean(ll), 4),
          stringsAsFactors = FALSE)
        pr_hd <- .prog_tick(pr_hd, sprintf("p=%d %s | possible=%d perfold=%d union=%.1f | FWERfalse=%.1f%% correct=%s",
          p, sprintf("n=%d %s",nn,HD_LABEL[[dg]]), p*(p-1L)/2L, choose(tf,2L), .safe_mean(nunion),
          out[[k]]$FWER_false_selection_pct,
          if(is.na(out[[k]]$Correct_pair_pct)) "n/a" else sprintf("%.1f%%", out[[k]]$Correct_pair_pct)))
      }}}
  .prog_done(pr_hd)
  do.call(rbind, out)}

sensitivity_analysis <- function(num_sim=100,n=500,base_seed=BASE_SEED){
  is_pair34<-function(f){if(f$n_inter<1)return(FALSE)
    any(vapply(f$inter,function(it)all(sort(c(it$j,it$k))==c(3,4)),TRUE))}
  out<-list();k<-0
  .n_blocks <- length(R2_SENS)*3L + length(C_SENS) +
               length(unique(c(100,n)))*length(M_SENS)*2L +
               length(LAMINT_SENS)*2L +
               4L*length(MININT_SENS)*2L +
               2L*3L
  pr_s <- .prog_new(.n_blocks, "Sensitivity analysis")
  .tick <- function() pr_s <<- .prog_tick(pr_s, sprintf("%s = %s",
              out[[k]]$Parameter, if(is.na(out[[k]]$Value)) "-" else as.character(out[[k]]$Value)))
  for(th in R2_SENS) for(dg in c("null_mixed","logit_interact","logit_pure_int")){
    anysel<-0;pair34<-0;nsucc<-0L;br<-c();lo<-c();au<-c()
    for(r in 1:num_sim){set.seed(base_seed+r);d<-generate_dataset(dg,n,0.25);sp<-split_stratified(d$y)
      f<-tryCatch(brcm_fit(d$X[sp$train,,drop=FALSE],d$y[sp$train],CAT_MASK,r2_thresh=th),error=function(e)NULL)
      if(is.null(f)||!isTRUE(f$converged))next                 
      nsucc<-nsucc+1L
      anysel<-anysel+(f$n_inter>0);pair34<-pair34+is_pair34(f)
      p<-brcm_predict(f,d$X[sp$test,,drop=FALSE]);yt<-d$y[sp$test]
      br<-c(br,mean((p-yt)^2));lo<-c(lo,.ll(yt,p));au<-c(au,.auc(yt,p))}
    k<-k+1;out[[k]]<-data.frame(Parameter="r2_thresh",Value=th,DGP=DGP_LABEL[dg],
      Attempted=num_sim,Successful=nsucc,
      Pct_any_interaction=if(nsucc)round(100*anysel/nsucc,1) else NA_real_,
      Pct_correct_pair=if(nsucc)round(100*pair34/nsucc,1) else NA_real_,
      Mean_AUC=round(.safe_mean(au),4),Mean_Brier=round(.safe_mean(br),4),
      Mean_LogLoss=round(.safe_mean(lo),4),Median_lambda=NA_real_,
      Mean_n_inter=NA_real_,Mean_ICI=NA_real_,stringsAsFactors=FALSE);.tick()}
  for(cs in C_SENS){
    nsucc<-0L;br<-c();lo<-c();au<-c();lam<-c()
    for(r in 1:num_sim){set.seed(base_seed+r);d<-generate_dataset("logit_mixed",n,0.25);sp<-split_stratified(d$y)
      f<-tryCatch(brcm_fit(d$X[sp$train,,drop=FALSE],d$y[sp$train],CAT_MASK,c_spear=cs),error=function(e)NULL)
      if(is.null(f)||!isTRUE(f$converged))next
      nsucc<-nsucc+1L
      p<-brcm_predict(f,d$X[sp$test,,drop=FALSE]);yt<-d$y[sp$test]
      br<-c(br,mean((p-yt)^2));lo<-c(lo,.ll(yt,p));au<-c(au,.auc(yt,p));lam<-c(lam,f$lambda)}
    k<-k+1;out[[k]]<-data.frame(Parameter="c_spearman",Value=cs,DGP="Logit (mixed)",
      Attempted=num_sim,Successful=nsucc,
      Pct_any_interaction=NA_real_,Pct_correct_pair=NA_real_,
      Mean_AUC=round(.safe_mean(au),4),Mean_Brier=round(.safe_mean(br),4),
      Mean_LogLoss=round(.safe_mean(lo),4),Median_lambda=round(.safe_median(lam),3),
      Mean_n_inter=NA_real_,Mean_ICI=NA_real_,stringsAsFactors=FALSE);.tick()}
  for(nn_m in unique(c(100, n))) for(mm in M_SENS) for(dg in c("logit_mixed","logit_nl")){
    nsucc<-0L;br<-c();lo<-c();au<-c();ec<-c()
    for(r in 1:num_sim){set.seed(base_seed+r);d<-generate_dataset(dg,nn_m,0.25);sp<-split_stratified(d$y)
      f<-tryCatch(brcm_fit(d$X[sp$train,,drop=FALSE],d$y[sp$train],CAT_MASK,m=mm),error=function(e)NULL)
      if(is.null(f)||!isTRUE(f$converged))next
      nsucc<-nsucc+1L
      pp<-brcm_predict(f,d$X[sp$test,,drop=FALSE]);yt<-d$y[sp$test]
      br<-c(br,mean((pp-yt)^2));lo<-c(lo,.ll(yt,pp));au<-c(au,.auc(yt,pp));ec<-c(ec,.ici(yt,pp))}
    k<-k+1;out[[k]]<-data.frame(Parameter="m_bernstein",Value=mm,
      DGP=paste0(DGP_LABEL[dg]," [n=",nn_m,"]"),
      Attempted=num_sim,Successful=nsucc,Pct_any_interaction=NA_real_,Pct_correct_pair=NA_real_,
      Mean_AUC=round(.safe_mean(au),4),Mean_Brier=round(.safe_mean(br),4),
      Mean_LogLoss=round(.safe_mean(lo),4),Median_lambda=NA_real_,
      Mean_n_inter=NA_real_,Mean_ICI=round(.safe_mean(ec),4),
      stringsAsFactors=FALSE);.tick()}
  for(li in LAMINT_SENS) for(dg in c("logit_interact","logit_pure_int")){
    nsucc<-0L;br<-c();lo<-c();au<-c();ni<-c()
    for(r in 1:num_sim){set.seed(base_seed+r);d<-generate_dataset(dg,n,0.25);sp<-split_stratified(d$y)
      f<-tryCatch(brcm_fit(d$X[sp$train,,drop=FALSE],d$y[sp$train],CAT_MASK,lam_int=li),error=function(e)NULL)
      if(is.null(f)||!isTRUE(f$converged))next
      nsucc<-nsucc+1L;ni<-c(ni,f$n_inter)
      pp<-brcm_predict(f,d$X[sp$test,,drop=FALSE]);yt<-d$y[sp$test]
      br<-c(br,mean((pp-yt)^2));lo<-c(lo,.ll(yt,pp));au<-c(au,.auc(yt,pp))}
    k<-k+1;out[[k]]<-data.frame(Parameter="lam_int",Value=li,DGP=DGP_LABEL[dg],
      Attempted=num_sim,Successful=nsucc,
      Pct_any_interaction=if(nsucc)round(100*mean(ni>0),1) else NA_real_,
      Pct_correct_pair=NA_real_,
      Mean_AUC=round(.safe_mean(au),4),Mean_Brier=round(.safe_mean(br),4),
      Mean_LogLoss=round(.safe_mean(lo),4),Median_lambda=NA_real_,
      Mean_n_inter=round(.safe_mean(ni),3),Mean_ICI=NA_real_,stringsAsFactors=FALSE);.tick()}
  for(nn_g in c(200,400,700,1200)) for(mi in MININT_SENS) for(dg in c("null_mixed","logit_interact")){
    nsucc<-0L;anysel<-0L;correct<-0L;au<-c();lo<-c()
    for(r in 1:num_sim){set.seed(base_seed+r);d<-generate_dataset(dg,nn_g,0.25);sp<-split_stratified(d$y)
      f<-tryCatch(brcm_fit(d$X[sp$train,,drop=FALSE],d$y[sp$train],CAT_MASK,min_n_inter=mi),error=function(e)NULL)
      if(is.null(f)||!isTRUE(f$converged))next
      nsucc<-nsucc+1L;anysel<-anysel+(f$n_inter>0)
      correct<-correct+(f$n_inter>0 && any(vapply(f$inter,function(it)all(sort(c(it$j,it$k))==c(3,4)),TRUE)))
      pp<-brcm_predict(f,d$X[sp$test,,drop=FALSE]);yt<-d$y[sp$test]
      au<-c(au,.auc(yt,pp));lo<-c(lo,.ll(yt,pp))}
    k<-k+1;out[[k]]<-data.frame(Parameter="min_n_inter",Value=mi,
      DGP=paste0(DGP_LABEL[dg]," [n=",nn_g,"]"),
      Attempted=num_sim,Successful=nsucc,
      Pct_any_interaction=if(nsucc)round(100*anysel/nsucc,1) else NA_real_,
      Pct_correct_pair=if(nsucc&&dg!="null_mixed")round(100*correct/nsucc,1) else NA_real_,
      Mean_AUC=round(.safe_mean(au),4),Mean_Brier=NA_real_,
      Mean_LogLoss=round(.safe_mean(lo),4),Median_lambda=NA_real_,
      Mean_n_inter=NA_real_,Mean_ICI=NA_real_,stringsAsFactors=FALSE);.tick()}
  for(sc in c("r2","logloss")) for(dg in c("null_mixed","logit_interact","logit_pure_int")){
    nsucc<-0L;anysel<-0L;correct<-0L;au<-c();lo<-c()
    for(r in 1:num_sim){set.seed(base_seed+r);d<-generate_dataset(dg,n,0.25);sp<-split_stratified(d$y)
      f<-tryCatch(brcm_fit(d$X[sp$train,,drop=FALSE],d$y[sp$train],CAT_MASK,sel_crit=sc),error=function(e)NULL)
      if(is.null(f)||!isTRUE(f$converged))next
      nsucc<-nsucc+1L;anysel<-anysel+(f$n_inter>0)
      correct<-correct+(f$n_inter>0 && any(vapply(f$inter,function(it)all(sort(c(it$j,it$k))==c(3,4)),TRUE)))
      pp<-brcm_predict(f,d$X[sp$test,,drop=FALSE]);yt<-d$y[sp$test]
      au<-c(au,.auc(yt,pp));lo<-c(lo,.ll(yt,pp))}
    k<-k+1;out[[k]]<-data.frame(Parameter="sel_crit",Value=NA_real_,DGP=paste0(DGP_LABEL[dg]," [",sc,"]"),
      Attempted=num_sim,Successful=nsucc,
      Pct_any_interaction=if(nsucc)round(100*anysel/nsucc,1) else NA_real_,
      Pct_correct_pair=if(nsucc&&dg!="null_mixed")round(100*correct/nsucc,1) else NA_real_,
      Mean_AUC=round(.safe_mean(au),4),Mean_Brier=NA_real_,
      Mean_LogLoss=round(.safe_mean(lo),4),Median_lambda=NA_real_,
      Mean_n_inter=NA_real_,Mean_ICI=NA_real_,stringsAsFactors=FALSE);.tick()}
  do.call(rbind,out)}

run_full_simulation <- function(num_sim=NUM_SIM,ns=AYAR_ns,prevs=AYAR_prev,
                                n_cores=N_CORES,out_dir=OUT_DIR,base_seed=BASE_SEED,
                                do_sensitivity=TRUE,sens_reps=NULL,
                                do_highdim=TRUE,hd_reps=NULL,hd_ps=c(10,20,50),hd_ns=c(500,1000),
                                hd_top_feats=8L){
  if(is.null(hd_reps)) hd_reps <- min(num_sim, 500L)
  if(is.null(sens_reps)) sens_reps <- min(num_sim, 200L)
  .check_pkgs();dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  registry<-build_registry()
  cat("Methods:",paste(names(registry),collapse=", "),"\n")
  cat("Calibrated variants:",paste(intersect(CALIBRATE,names(registry)),collapse=", "),"\n")
  cl<-NULL
  if(n_cores>1&&has("parallel")){
    cl<-parallel::makeCluster(n_cores)
    parallel::clusterExport(cl,varlist=ls(envir=.GlobalEnv,all.names=TRUE),envir=.GlobalEnv)
    parallel::clusterEvalQ(cl,{suppressMessages({
      for(p in c("mgcv","rpart","randomForest","e1071","xgboost","lightgbm","catboost"))
        if(requireNamespace(p,quietly=TRUE)) library(p,character.only=TRUE)});NULL})
    on.exit(try(parallel::stopCluster(cl),silent=TRUE),add=TRUE)}
  raw_all<-list();diag_all<-list();t0<-Sys.time()
  n_scen<-length(BASE_DGPS)*length(ns)*length(prevs)
  ckdir<-file.path(out_dir,"_checkpoints");dir.create(ckdir,recursive=TRUE,showWarnings=FALSE)
  ck_file<-function(scn) file.path(ckdir,paste0(scn,".rds"))
  mf<-file.path(ckdir,"_manifest.rds")
  cfg<-list(num_sim=num_sim, ns=ns, prevs=prevs, dgps=BASE_DGPS,
            base_seed=base_seed, n_methods=length(registry))
  if(file.exists(mf)){
    old_cfg<-tryCatch(readRDS(mf),error=function(e)NULL)
    if(!is.null(old_cfg)){
      diffs<-character(0)
      cmp<-function(nm,a,b) if(!identical(a,b))
        sprintf("%s: checkpoints were made with %s, this run uses %s",
                nm, paste(a,collapse=","), paste(b,collapse=","))
      for(k in names(cfg)){
        d<-cmp(k, old_cfg[[k]], cfg[[k]]); if(length(d)&&nzchar(d)) diffs<-c(diffs,d)}
      if(length(diffs))
        stop("The existing checkpoints in ", ckdir,
             " do not match this run:\n  - ", paste(diffs, collapse="\n  - "),
             "\nUse a different BRCM_OUT, or delete that directory to start fresh.")}
  } else saveRDS(cfg, mf)
  failed_scen<-character(0)
  done_before<-length(list.files(ckdir,pattern="\\.rds$"))
  if(done_before)
    cat(sprintf("Resuming: %d of %d scenarios already on disk in %s\n",
                done_before,n_scen,ckdir))
  pr<-.prog_new(n_scen,sprintf("Factorial study (%d reps/scenario)",num_sim))
  for(dg in BASE_DGPS) for(nn in ns) for(pv in prevs){
    scn<-sprintf("%s_n%d_p%02d",dg,nn,round(pv*100))
    f<-ck_file(scn)
    if(file.exists(f)){
      r<-tryCatch(readRDS(f),error=function(e)NULL)          
      if(!is.null(r)){
        if(isTRUE(r$failed)){failed_scen<-c(failed_scen,scn)
          pr<-.prog_tick(pr,sprintf("%s (previously failed, skipped)",scn));next}
        raw_all[[scn]]<-r$raw;diag_all[[scn]]<-r$diag
        pr<-.prog_tick(pr,sprintf("%s (from checkpoint)",scn));next}}
    r<-tryCatch(run_scenario(dg,nn,pv,num_sim,registry,base_seed,cl,
                             progress=TRUE,progress_label=sprintf("%s reps",scn),
                             n_cores=n_cores),
                error=function(e){cat(sprintf("\n  [!] %s failed: %s\n      rebuilding the cluster and retrying once\n",
                                              scn,conditionMessage(e)));NULL})
    if(!is.null(r) && !identical(r$cl,cl)) cl<-r$cl   
    if(is.null(r)){
      if(!is.null(cl)) try(parallel::stopCluster(cl),silent=TRUE)
      cl<-tryCatch({k<-parallel::makeCluster(n_cores)
        parallel::clusterExport(k,varlist=ls(envir=.GlobalEnv,all.names=TRUE),envir=.GlobalEnv)
        parallel::clusterEvalQ(k,{suppressMessages({library(stats);library(utils)})});k},
        error=function(e){cat("      cluster rebuild failed; continuing sequentially\n");NULL})
      r<-tryCatch(run_scenario(dg,nn,pv,num_sim,registry,base_seed,cl,
                               progress=TRUE,progress_label=sprintf("%s retry",scn),
                               n_cores=n_cores),
                  error=function(e){cat(sprintf("  [!] %s failed again: %s\n",scn,conditionMessage(e)));NULL})
      if(!is.null(r) && !identical(r$cl,cl)) cl<-r$cl}
    if(is.null(r)){
      failed_scen<-c(failed_scen,scn)
      saveRDS(list(raw=list(),diag=NULL,failed=TRUE),ck_file(scn))   
      pr<-.prog_tick(pr,sprintf("%s (SKIPPED after failure)",scn));next}
    r$diag$Scenario<-scn
    tmp<-paste0(f,".tmp");saveRDS(list(raw=r$raw,diag=r$diag),tmp);file.rename(tmp,f)
    raw_all[[scn]]<-r$raw;diag_all[[scn]]<-r$diag
    pr<-.prog_tick(pr,sprintf("%s (%d methods)",scn,length(r$raw)))}
  .prog_done(pr)
  if(length(failed_scen)){
    cat("\n",strrep("!",70),"\n",sep="")
    cat(sprintf(" %d scenario(s) could not be completed and were skipped:\n",length(failed_scen)))
    for(z in failed_scen) cat("   ",z,"\n")
    cat(" All other results are complete. To retry one of them after fixing the\n")
    cat(" cause, delete its file from ",ckdir," and re-run.\n",sep="")
    cat(strrep("!",70),"\n\n",sep="")}
  P<-paired_analysis(raw_all)
  idg<-NULL
  if(do_sensitivity) idg<-tryCatch(sensitivity_analysis(sens_reps),error=function(e){
    if(isTRUE(STRICT_PKGS)) stop("sensitivity_analysis failed in a strict run: ",conditionMessage(e))
    warning("sensitivity_analysis failed: ",conditionMessage(e));NULL})
  hd<-NULL
  if(do_highdim){
    cat("High-dimensional screening analysis ...\n")
    hd<-tryCatch(highdim_analysis(num_sim=hd_reps,ps=hd_ps,ns=hd_ns,top_feats=hd_top_feats),error=function(e){
      if(isTRUE(STRICT_PKGS)) stop("highdim_analysis failed in a strict run: ",conditionMessage(e))
      cat(" (failed:",conditionMessage(e),")\n");NULL})}
  sim<-list(raw=raw_all,paired=P,diag=do.call(rbind,diag_all),inter_diag=idg,highdim=hd,
    config=list(num_sim=num_sim,ns=ns,prevs=prevs,oof_K=OOF_K,r2_thresh=R2_THRESH,
      base_seed=base_seed,test_frac=TEST_FRAC,methods=names(registry),
      sessionInfo=utils::capture.output(utils::sessionInfo())))
  saveRDS(sim,file.path(out_dir,"sim.rds"))
  write_tables(sim,out_dir);write_figures(sim,out_dir)
  cat("Elapsed:",format(round(difftime(Sys.time(),t0),1)),"| output in",normalizePath(out_dir),"\n")
  invisible(sim)}

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

brcm_diagnose <- function(dgp, n = 100, prev = 0.50, reps = 300, registry = NULL){
  if(is.null(registry)) registry <- build_registry()
  cat(sprintf("Scenario %s, n = %d, prevalence = %.2f, %d replications\n",
              dgp, n, prev, reps))
  cat(sprintf("Testing %d methods, each alone in its own worker process\n\n", length(registry)))
  out <- data.frame(Method = names(registry), Status = NA_character_,
                    Successful = NA_integer_, stringsAsFactors = FALSE)
  for(i in seq_along(registry)){
    nm <- names(registry)[i]
    cat(sprintf("  %-18s ... ", nm)); flush.console()
    cl <- parallel::makeCluster(1L)
    parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv, all.names = TRUE),
                            envir = .GlobalEnv)
    r <- tryCatch({
      k <- parallel::parLapply(cl, 1L, function(z){
        one <- registry[[nm]]; n_ok <- 0L
        for(rep_id in seq_len(reps)){
          set.seed(BASE_SEED + rep_id)
          d  <- generate_dataset(dgp, n, prev)
          cm <- dgp_cat_mask(dgp)
          sp <- split_stratified(d$y)
          pp <- try(one(d$X[sp$train,,drop=FALSE], d$y[sp$train],
                        d$X[sp$test,,drop=FALSE], cm), silent = TRUE)
          if(!inherits(pp,"try-error") && all(is.finite(pp))) n_ok <- n_ok + 1L}
        n_ok})[[1]]
      list(st = "ok", n = k)
    }, error = function(e) list(st = paste("CRASHED:", conditionMessage(e)), n = NA_integer_))
    try(parallel::stopCluster(cl), silent = TRUE)
    out$Status[i] <- r$st; out$Successful[i] <- r$n
    cat(r$st, if(!is.na(r$n)) sprintf(" (%d/%d fits)", r$n, reps) else "", "\n", sep = "")}
  bad <- out$Method[grepl("^CRASHED", out$Status)]
  cat("\n", strrep("=", 60), "\n", sep = "")
  print(out, row.names = FALSE)
  cat(strrep("=", 60), "\n", sep = "")
  if(length(bad)) cat("\n>>> CULPRIT:", paste(bad, collapse = ", "), "\n")
  else            cat("\n>>> No method crashed alone; try a larger `reps`.\n")
  invisible(out)}

brcm_run_all <- function(){
  if(dir.exists(OUT_DIR) && !dir.exists(file.path(OUT_DIR,"_checkpoints"))){
    kept <- paste0(OUT_DIR, "_previous_", format(Sys.Date(), "%Y%m%d"))
    if(file.rename(OUT_DIR, kept))
      cat("Existing output directory preserved as:", kept, "\n")}

  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("\nOutput directory: %s   (%d replications per scenario)\n",
              normalizePath(OUT_DIR), NUM_SIM))
  cat("\nSelf-test\n")
  if(!brcm_self_test(verbose = TRUE))
    stop("Self-test failed. The run was not started; the code is not in a usable state.")

  cat(sprintf("\nFactorial study: %d scenarios x %d replications\n",
              length(BASE_DGPS)*length(AYAR_ns)*length(AYAR_prev), NUM_SIM))
  sim <- run_full_simulation(do_highdim = FALSE, do_sensitivity = FALSE)
  saveRDS(sim, file.path(OUT_DIR, "01_factorial.rds"))
  cat(">>> factorial complete, saved to", file.path(OUT_DIR, "01_factorial.rds"), "\n")

  cat("\nHigh-dimensional screening analysis\n")
  hd <- highdim_analysis(num_sim = min(NUM_SIM, 500L), ps = c(10, 20, 50),
                         ns = c(500, 1000), top_feats = 8L)
  saveRDS(hd, file.path(OUT_DIR, "02_highdim.rds"))
  cat(">>> high-dimensional analysis complete\n")

  cat("\nSensitivity analysis\n")
  sens <- sensitivity_analysis(num_sim = min(NUM_SIM, 200L))
  saveRDS(sens, file.path(OUT_DIR, "03_sensitivity.rds"))
  cat(">>> sensitivity analysis complete\n")

  sim$highdim <- hd; sim$inter_diag <- sens
  saveRDS(sim, file.path(OUT_DIR, "sim.rds"))
  write_tables(sim, OUT_DIR); write_figures(sim, OUT_DIR)
  td <- file.path(OUT_DIR, "tables"); dir.create(td, recursive = TRUE, showWarnings = FALSE)
  if(!is.null(hd))   utils::write.csv(hd,   file.path(td, "T_highdim_screening.csv"), row.names = FALSE)
  if(!is.null(sens)) utils::write.csv(sens, file.path(td, "T_interaction_diagnostics.csv"), row.names = FALSE)
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("  COMPLETE\n")
  cat("  directory :", normalizePath(OUT_DIR), "\n")
  cat("  tables    :", length(list.files(file.path(OUT_DIR,"tables"))), "files\n")
  cat("  figures   :", length(list.files(file.path(OUT_DIR,"figures"))), "files\n")
  cat("  scenarios :", length(sim$raw), "\n")
  cat(strrep("=", 60), "\n", sep = "")
  invisible(sim)}

if(AUTORUN) sim <- brcm_run_all()
