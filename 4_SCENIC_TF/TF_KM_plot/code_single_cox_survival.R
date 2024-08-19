library(survival)
library(survminer)
df = read.table("./TF_ESCA_survival_data.txt", row.names=1, header=T, sep="\t")
P_arr = c()
P_arr[1] = 0
P_arr[2] = 0
for (i in 3:ncol(df)){
	sub_df = df[, c(1, 2, i)]
	sub_df = na.omit(sub_df)
	names(sub_df)<-c("time", "status", "Site1")
	result = coxph(Surv(time, status)~Site1 , data=sub_df)
	p_value = summary(result)$coef[1,5]
	P_arr[i] = p_value
}
P_df = data.frame(P=P_arr)
rownames(P_df) = colnames(df)
write.table(P_df, "./single_cox_p_value.txt", sep="\t", quote=F)

sub_df = subset(df, select=c("time", "status", "DDIT3"))
median = median(sub_df$DDIT3)
sub_df["DDIT3_Type"] = ifelse(sub_df$DDIT3>median, 1, 0)

sub_df = read.table("./sub_df.txt", header=T, row.names=1, sep="\t")
fit = survfit(Surv(time, status)~DDIT3_Type, data=sub_df)
ggsurvplot(fit,
	data=sub_df,
	conf.int=TRUE,
	pval=TRUE,
	surv.median.line="hv",
	risk.table=TRUE,
	xlab="Time(days)",
	palette=c("#0da9ce", "#e74a32"),
	legend.title="",
	legend.labs=c("Low expression", "High expression"),
	break.x.by=500,
	pval.size=7,
	font.legend=13
)


