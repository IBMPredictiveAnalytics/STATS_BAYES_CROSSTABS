#/***********************************************************************
# * Licensed Materials - Property of IBM 
# *
# * IBM SPSS Products: Statistics Common
# *
# * (C) Copyright IBM Corp. 2015
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp. 
# ************************************************************************/

# author__ = "IBM SPSS, JKP"
# version__ = "1.0.0"

# History
# 05-oct-2015 Original Version


gtxt <- function(...) {
    return(gettext(...,domain="STATS_BAYES_CROSSTABS"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_BAYES_CROSSTABS"))
}

sampletypes=list(poisson="poisson", jointmulti="jointMulti",
    indepmulti="indepMulti", hypergeometric="hypergeom")
fixedmargins=list(rows="rows", columns="cols")

### MAIN ROUTINE ###
doBayesxtab = function(rows, columns, sampletype="poisson", fixedmargin="rows",
    showpost=TRUE, priorconc=1., iterations=1000,
    displaytable=TRUE, workspaceaction="clear", modelfileout=NULL) {
    # Estimate Bayes crosstabulation
    
    # The modelsource and modelfile
    # parameters are not implemented, awaiting requests for that functionality

    setuplocalization("STATS_BAYES_CROSSTABS")
    
    # A warnings proc name is associated with the regular output
    # (and the same omsid), because warnings/errors may appear in
    # a separate procedure block following the regular output
    procname=gtxt("Bayesian Crosstab")
    warningsprocname = gtxt("Bayesian Crosstab: Warnings")
    omsid="STATSBAYESCROSSTAB"
    warns = Warn(procname=warningsprocname,omsid=omsid)

    tryCatch(library(BayesFactor), error=function(e){
        warns$warn(gtxtf("The R %s package is required but could not be loaded.", "BayesFactor"),dostop=TRUE)
        }
    )
    spsswtvar = spssdictionary.GetWeightVariable()

    if (!is.null(spssdata.GetSplitVariableNames())) {
        warns$warn(
            gtxt("Split variables are not honored by this procedure"),
            dostop=FALSE)
    }
    if (sampletype == "hypergeometric" && showpost) {
        warns$warn(gtxt("Posterior distribution statistics are not available for hypergeometric sampling."),
            dostop=FALSE)
        showpost = FALSE
    }

    alldata = c(rows, columns, spsswtvar)
    allargs = as.list(environment())
    # with only one variable, value is not a data frame
    dta = spssdata.GetDataFromSPSS(alldata, missingValueToNA=TRUE,
        factorMode="labels")
    if (!(is.factor(dta[[rows]]) && is.factor(dta[[columns]]))) {
        warns$warn(gtxt("The row and column variables must be categorical"),
            dostop=TRUE)
    }
    if (!is.null(spsswtvar)) {
        if (is.factor(dta[[spsswtvar]])) {
            warns$warn(gtxt("The SPSS weight cannot be categorical"),
                dostop=TRUE)
        }
    }

    # The procedure does not allow missing values
    allargs$ncases = nrow(dta)
    # If the data frame has only one column, the extraction of
    # complete cases loses the data frame class :-(
    if (length(dta) == 1) {
        dta[2] = 1:allargs$ncases
    }
    dta = dta[complete.cases(dta),]
    allargs$nvalid = nrow(dta)
    if (!is.null(spsswtvar)) {
        allargs$wtcount = sum(dta[[spsswtvar]])
    } else {
        allargs$wtcount = allargs$nvalid
    }
    # accumulate crosstab table and, if there is an SPSS weight
    # 

    frml = buildformula(rows, columns, spsswtvar)
    tt = xtabs(frml, data=dta, na.action=na.omit)
    if (sampletype == "hypergeometric" && (nrow(tt) != 2 || ncol(xt) != 2)) {
        warns$warn(gtxt("The hypergeometric sample type can only be used with 2 x 2 tables"),
            dostop=TRUE)
    }
    if (!is.null(spsswtvar)) {
        for (i in length(tt)) {
            tt[[i]] = round(tt[[i]])
        }
    }
    allargs$tt = tt
    res = tryCatch(
        contingencyTableBF(tt, sampleType=sampletypes[[sampletype]], 
            fixedMargin=fixedmargins[[fixedmargin]], priorConcentration = priorconc),
        error=function(e) {
            warns$warn(e, dostop=TRUE)
        }
    )
    post = doposterior(allargs, res, warns)

    displayresults(allargs, res, post, warns)
    
    if (!is.null(modelfileout)) {
        save(allargs, res, post, file=modelfileout)
    }
    if (workspaceaction == "retain") {
        assign("allargs", allargs, envir=.GlobalEnv)
        assign("res", res, envir=.GlobalEnv)
        assign("post", post, envir=.GlobalEnv)
    }
    warns$display()
}

buildformula = function(rows, columns, spsswtvar) {
    # return a forumla object
    if (is.null(spsswtvar)) {
        return(sprintf("~%s+%s", rows, columns))
    } else {
        return(sprintf("%s~%s+%s", spsswtvar, rows, columns))
    }
}

doposterior = function(allargs, res, warns) {
    # calculate posterior distribution if model index specified
    
    if (!allargs$showpost) {
        return(NULL)
    }

    post = posterior(model=res,iterations=allargs$iterations,
        progress=FALSE)
    return(post)
}

waction=list("clear"="clear", "retain"="retain")
outputmap = list("delta"=gtxt("Effect Size"), 
                 "mu"=gtxt("Mean"),
                 "beta (x - y)"=gtxt("Difference"),
                 "sig2"=gtxt("Sigma2"),
                 "g" = gtxt("g")
)
sampletypesTranslation = list(poisson=gtxt("Poisson"), jointmulti=gtxt("Joint Multinomial"),
    indepmulti=gtxt("Independent Multinomial"), hypergeometric=gtxt("Hypergeometric")
)
fixedmarginsTranslation=list(rows=gtxt("rows"), columns=gtxt("columns")
)

displayresults = function(allargs, res, post, warns) {
    # display results
    # allargs is the parameter set

    StartProcedure(allargs[["procname"]], allargs[["omsid"]])
    
    # summary results
    # input specifications
    # although groups can be specified (cengroup), separate results are not
    # produced.
    lbls = c(gtxt("Row Variable"),
             gtxt("Column Variable"),
             gtxt("Sample Type"),
             gtxt("Fixed Margin"),
             gtxt("Prior Concentration"),
             gtxt("Number of Cases"),
             gtxt("Number of Valid Cases"),
             gtxt("Weight Variable"),
             gtxt("Weighted Number of Valid Cases"),
             gtxt("Posterior Iterations"),
             gtxt("Workspace Action"),
             gtxt("Output Model File")
    )
    vals = c(
            allargs$rows,
            allargs$columns,
            sampletypesTranslation[[allargs$sampletype]],
            ifelse(allargs$sampletype == "indepmulti", 
                fixedmarginsTranslation[[allargs$fixedmargin]], gtxt("--NA--")),
            allargs$priorconc,
            allargs$ncases,
            allargs$nvalid,
            ifelse(is.null(allargs$spsswtvar), gtxt("--NA--"), allargs$spsswtvar),
            round(allargs$wtcount, 3),
            ifelse(!allargs$showpost, gtxt("--NA--"), allargs$iterations),
            waction[allargs$workspaceaction],
            ifelse(is.null(allargs$modelfile), gtxt("--NA--"), allargs$modelfile)
    )
    

    spsspivottable.Display(data.frame(cbind(vals), row.names=lbls), title = gtxt("Summary"),
        collabels=c(gtxt("Summary")), templateName="BAYESCROSSTABSSUMMARY", outline=gtxt("Bayes Crosstabs Summary"),
        caption = gtxtf("Computations done by R package BayesFactor, version: %s", packageVersion("BayesFactor"))
    )
    if (allargs$displaytable) {
        table = as.data.frame.matrix(allargs$tt)
        # Add marginal totals.  Leading blank in total labels is to ensure uniqueness of labels
        table = cbind(table, apply(table, 1, sum, na.rm=TRUE))
        table = rbind(table, apply(table, 2, sum, na.rm=TRUE))
        names(table)[ncol(table)] = gtxt(" Total")
        row.names(table)[nrow(table)] = gtxt(" Total")
        
        spsspivottable.Display(
            table,
            templateName="BAYESCROSSTABSTABLE", rowdim=allargs$rows, coldim=allargs$columns,
            hiderowdimtitle=FALSE, hidecoldimtitle=FALSE,
            title=gtxtf("Crosstabulation: %s by %s", allargs$rows, allargs$columns),
            outline=gtxt("Crosstabulation")
        )
    }

    ressum = extractBF(res)
    bf = data.frame(seq(1: length(res)),ressum[1:2])
    if (is.nan(bf$bf)) {
        warns$warn(gtxt("The Bayes factor could not be computed.  You may need to choose a different prior concentration value."),
            dostop=TRUE)
    }
    bf[3] = bf[3] * 100.
    bf = data.frame(bf, length(res) - rank(bf[2]) + 1)
    # add in posterior probabilities excluding Intercept only model
    ###postprob = data.frame(as.BFprobability(newPriorOdds(res) * res))[-(nrow(bf)+1), 1]
    
    # construct posterior probabilities and merge with Bayes factors
    # The order for probabilities may not be the same as for the Bayes factors
    # which requires a few extra steps to get things merged
    # the BF data frame may not have the intercept row, so that row may be discarded
    postprob = data.frame(as.BFprobability(newPriorOdds(res) * res))[1]

    bf = merge(bf, postprob, by="row.names")
    bf = bf[order(bf[[2]]),]
    row.names(bf) = bf[["Row.names"]]
    bf = bf[-1]
    
    bf = bf[c(2,3,5)]

    names(bf) = c(
        gtxt("Bayes Factor"), gtxt("Error (+-%)"),
        gtxt("Posterior Probabilities (Equal Prior)"))
    row.names(bf) = c("Non-independence")
    spsspivottable.Display(bf,
        title=gtxt("Bayes Factors"),
        rowdim=gtxt("Prior Scale"), 
        hiderowdimtitle=FALSE,
        templateName="BAYESCROSSTABSFACTORS",
        outline=gtxt("Bayes Factors")
    )

    if (allargs$showpost) {
        # construct row names to match the xtab table labels
        labels=list()
        i = 1
        dd = dimnames(allargs$tt)
        dd1 = dd[[1]]
        dd2 = dd[[2]]
        for (x2 in dd2) {
            for (x1 in dd1) {
                labels[i] = paste(x1, x2, sep=", ")
                i = i+1
            }
        }

        postsum = summary(post)
        postsumstats = data.frame(postsum$statistics[,-4])  # omit time series SEs
        row.names(postsumstats) = labels
        names(postsumstats) = c(gtxt("Mean"), gtxt("Std. Deviation"), gtxt("SE Mean"))
        spsspivottable.Display(
            postsumstats, 
            title=gtxt("Posterior Summary Statistics"),
            rowdim=gtxt("Statistic"),
            hiderowdimtitle=TRUE,
            templateName="BAYESCROSSTABSPOSTSTATS",
            outline=gtxt("Posterior Summary Statistics")
        )
        
        postsumquant = postsum$quantiles
        row.names(postsumquant) = labels
        spsspivottable.Display(
            postsumquant,
            title=gtxtf("Posterior Quantiles"),
            rowdim=gtxt("Variables"),
            hiderowdimtitle=TRUE,
            templateName="BAYESTTESTPOSTQUANTILES",
            outline=gtxt("Posterior Quantiles")
        )
    }

    spsspkg.EndProcedure()
}

Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = mylist2env(lcl) # makes this list into an environment

    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.

        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 

        if (is.null(msg) || dostop) {
            lcl$display(inproc)  # display messages and end procedure state
            if (dostop) {
                stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
            }
        }
    }
    
    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any

        if (lcl$msgnum == 0) {   # nothing to display
            if (inproc) {
                spsspkg.EndProcedure()
            }
        } else {
            if (!inproc) {
                procok =tryCatch({
                    StartProcedure(lcl$procname, lcl$omsid)
                    TRUE
                    },
                    error = function(e) {
                        FALSE
                    }
                )
            }
            if (procok) {  # build and display a Warnings table if we can
                table = spss.BasePivotTable("Warnings ","Warnings") # do not translate this
                rowdim = BasePivotTable.Append(table,Dimension.Place.row, 
                    gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)

                for (i in 1:lcl$msgnum) {
                    rowcategory = spss.CellText.String(as.character(i))
                    BasePivotTable.SetCategories(table,rowdim,rowcategory)
                    BasePivotTable.SetCellValue(table,rowcategory, 
                        spss.CellText.String(lcl$msglist[[i]]))
                }
                spsspkg.EndProcedure()   # implies display
            } else { # can't produce a table
                for (i in 1:lcl$msgnum) {
                    print(lcl$msglist[[i]])
                }
            }
        }
    }
    return(lcl)
}

mylist2env = function(alist) {
    env = new.env()
    lnames = names(alist)
    for (i in 1:length(alist)) {
        assign(lnames[[i]],value = alist[[i]], envir=env)
    }
    return(env)
}

# localization initialization
setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 
# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
        spsspkg.StartProcedure(procname, omsid)
    }
    else {
        spsspkg.StartProcedure(omsid)
    }
}



Run = function(args) {
    #Execute the STATS BAYES CROSSTAB command

    cmdname = args[[1]]
    args = args[[2]]
    oobj = spsspkg.Syntax(list(
        spsspkg.Template("ROWS", subc="", ktype="existingvarlist", var="rows"),
        spsspkg.Template("COLUMNS", subc="", ktype="existingvarlist", var="columns"),
        spsspkg.Template("SAMPLETYPE", subc="", ktype="str", var="sampletype",
            vallist=list("poisson", "jointmulti", "indepmulti", "hypergeometric")),
        spsspkg.Template("FIXEDMARGIN", subc="", ktype="str", var="fixedmargin",
            vallist=list("rows", "columns")),

        spsspkg.Template("POSTERIOR", subc="OPTIONS", ktype="bool", var="showpost"),
        spsspkg.Template("ITERATIONS", subc="OPTIONS", ktype="int", var="iterations",
            vallist=list(2)),
        spsspkg.Template("PRIORCONCENTRATION", subc="OPTIONS", ktype="float", var="priorconc",
            vallist=list(.0001)),
        spsspkg.Template("DISPLAYTABLE", subc="OPTIONS", ktype="bool", var="displaytable"),
        
        spsspkg.Template("WORKSPACE", subc="SAVE", ktype="str", var="workspaceaction",
            vallist=list("retain", "clear")),
        spsspkg.Template("MODELFILE", subc="SAVE", ktype="literal", var="modelfileout")
    ))

    # A HELP subcommand overrides all else
    if ("HELP" %in% attr(args,"names")) {
        helper(cmdname)
    }
    else {
        res <- spsspkg.processcmd(oobj, args, "doBayesxtab")
    }
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}
