calculate_macroF1 <- function(true_celltypes, pred_celltypes) {
    ## calculate macroF1
    union_celltypes = union(true_celltypes, pred_celltypes)
    cm = matrix(0, nrow=length(union_celltypes), ncol=length(union_celltypes))
    rownames(cm) = union_celltypes
    colnames(cm) = union_celltypes

    ## correct levels of target object
    predict_tb = table(true_celltypes, pred_celltypes)
    cm[rownames(predict_tb), colnames(predict_tb)] = predict_tb

    diag = diag(cm)
    precision = diag / colSums(cm) 
    precision[is.nan(precision)] = 0
    recall = diag / rowSums(cm) 
    recall[is.nan(recall)] = 0
    f1 = 2 * precision * recall / (precision + recall) 
    f1[is.nan(f1)] = 0
    macroF1 = mean(f1)
    return(macroF1)
}

calculate_accuracy <- function(true_celltypes, pred_celltypes) {
    return(sum(true_celltypes == pred_celltypes) / length(true_celltypes))
}

calculate_ARI <- function(true_celltypes, pred_celltypes) {
    require(aricode)
    return(ARI(true_celltypes, pred_celltypes))
}

