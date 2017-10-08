library(metabaser)
metabase.connect(dbname = "prd380", uid = "mb_user", pwd = "Temp_usr_9999", host = "metabase.am.lilly.com", port = 1525, type = "sid", driver = "jdbc")



drug_list <-
function(gene)
{
	gene_id = recognize.nwobj.names(gene)
	down <- get.nwobj.noodletable(start = gene_id)
	up <- get.nwobj.noodletable(end = gene_id)
	gene_down <- do.call("rbind", lapply(unique(down$noodle_id), get.noodle.genes))$gene
	gene_up <- do.call("rbind", lapply(unique(up$noodle_id), get.noodle.genes))$gene
	genes <- unique(c(gene_down, gene_up))
	drugs <- get.gene.drugs(genes)
	result <- data.frame(
		DrugName = drugs$drug_name,
		DrugTarget = drugs$target_name,
		DiseaseName = drugs$disease_name,
		Status = drugs$status
		)
	return(result)
}


gene_list <-
function(gene)
{
	gene_id = recognize.nwobj.names(gene)
	down <- get.nwobj.noodletable(start = gene_id)
	up <- get.nwobj.noodletable(end = gene_id)
	gene_down <- do.call("rbind", lapply(unique(down$noodle_id), get.noodle.genes))$genesymbol
	gene_up <- do.call("rbind", lapply(unique(up$noodle_id), get.noodle.genes))$genesymbol
	genes <- unique(c(gene_down, gene_up))
	drugs = character()
	for(i in 1: length(genes))
	{
		g <- recognize.gene.symbols(genes[i])
		drugs <- c(drugs, paste(unique(get.gene.drugs(g)$drug_name), collapse = ","))
	}	
	result <- data.frame(GeneSymbol = genes, Drugs = drugs)
	return(result)
}

noodles_detail <-
function(gene)
{
	#create an empty frame that will hold the result
	result <- data.frame(
		gene1 = character(),
		gene2 = character(),
		direction = character(),
		noodle = character(),
		effect = character(),
		chain = character(),
		distance = numeric(),
		effect_chain = character(),
		maps = character(),
		genes_in_get_noodle_genes = character(),
		genes_number = numeric()
		)
	gene_id = recognize.nwobj.names(gene)
	
	down <- get.nwobj.noodletable(start = gene_id)
	up <- get.nwobj.noodletable(end = gene_id)
	   
	#process down sequence
	if(nrow(down) > 0)
	{
		dis <- 0
		seq <- ""
		eff <- ""
		map <- ""
		genes_in_get_noodle_genes <- ""
		genes_number <- 0
		chn <- character()
		chne <- character()
		g_in_n <- ""
		g_n <- 0
		id <- ""  
		for(i in 1:length(down$noodle_id))
		{
			if(down$noodle_id[i] != id)
			{
				if(i > 1)
				{
					fr_down <- data.frame(
						gene1 = gene,
						gene2 = chn[length(chn)],
						direction = "down",
						noodle = seq,
						effect = eff,
						chain = paste(chn, collapse="--"),
						distance = dis - 1,
						effect_chain = paste(chne, collapse="|"),
						maps = map,
						genes_in_get_noodle_genes = g_in_n,
						genes_number = g_n
					)
					result <- rbind(result, fr_down)
					chn <- character()
					chne <- character()
					dis <- 0
				}
				id <- down$noodle_id[i]
				seq <- describe.noodles(id)$noodle_name
				#dis <- nchar(seq) - nchar(gsub("--", "-", seq))
				if(has.noodle.genes(id) == TRUE)
				{
					genes <- get.noodle.genes(id)$genesymbol
					g_n <- length(genes)
					g_in_n <- paste(genes, collapse = '--')
				}
				else
				{
					g_n <- 0
					g_in_n <- ""
				}
				if(length(get.noodle.maps(id)$map_name) > 0)
					map <- paste(get.noodle.maps(id)$map_name, collapse = '--')
				else
					map <- ""
				eff <- noodle.effect(describe.noodles(id))
				chn <- c(chn, describe.nwobjs(down$id[i])$name)
				chne <- c(chne, down$effect[i])
				dis <- dis + 1
			}
			else
				chn <- c(chn, describe.nwobjs(down$id[i])$name)
				chne <- c(chne, down$effect[i])
				dis <- dis + 1
			if (i == length(down$noodle_id))
			{
				fr_down <- data.frame(
					gene1 = gene,
					gene2 = chn[length(chn)],
					direction = "down",
					noodle = seq,
					effect = eff,
					chain = paste(chn, collapse="--"),
					distance = dis -1,
					effect_chain = paste(chne, collapse="|"),
					maps = map,
					genes_in_get_noodle_genes = g_in_n,
					genes_number = g_n
				)
				result <- rbind(result, fr_down)
			}

		}
	}
	
	#process up sequence
	if(nrow(up) > 0)
	{
		dis <- 0
		seq <- ""
		eff <- ""
		map <- ""
		genes_in_get_noodle_genes <- ""
		genes_number <- 0
		chn <- character()
		chne <- character()
		g_in_n <- ""
		g_n <- 0
		id <- ""  
		for(i in 1:length(up$noodle_id))
		{
			if(up$noodle_id[i] != id)
			{
				if(i > 1)
				{
					fr_up <- data.frame(
						gene1 = chn[1],
						gene2 = gene,
						direction = "up",
						noodle = seq,
						effect = eff,
						chain = paste(chn, collapse="--"),
						distance = dis - 1,
						effect_chain = paste(chne, collapse="|"),
						maps = map,
						genes_in_get_noodle_genes = g_in_n,
						genes_number = g_n
					)
					result <- rbind(result, fr_up)
					chn <- character()
					chne <- character()
					dis <- 0
				}
				id <- up$noodle_id[i]
				seq <- describe.noodles(id)$noodle_name
				#dis <- nchar(seq) - nchar(gsub("--", "-", seq))
				if(has.noodle.genes(id) == TRUE)
				{
					genes <- get.noodle.genes(id)$genesymbol
					g_n <- length(genes)
					g_in_n <- paste(genes, collapse = '--')
				}
				else
				{
					g_n <- 0
					g_in_n <- ""
				}
				if(length(get.noodle.maps(id)$map_name) > 0)
					map <- paste(get.noodle.maps(id)$map_name, collapse = '--')
				else
					map <- ""
				eff <- noodle.effect(describe.noodles(id))
				chn <- c(chn, describe.nwobjs(up$id[i])$name)
				chne <- c(chne, up$effect[i])
				dis <- dis + 1
			}
			else
				chn <- c(chn, describe.nwobjs(up$id[i])$name)
				chne <- c(chne, up$effect[i])
				dis <- dis + 1
			if (i == length(up$noodle_id))
			{
				fr_up <- data.frame(
					gene1 = chn[1],
					gene2 = gene,
					direction = "up",					
					noodle = seq,
					effect = eff,
					chain = paste(chn, collapse="--"),
					distance = dis -1 ,
					effect_chain = paste(chne, collapse="|"),
					maps = map,
					genes_in_get_noodle_genes = g_in_n,
					genes_number = g_n
				)
				result <- rbind(result, fr_up)
			}
		}
	}
return(result)
}

