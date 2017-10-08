library(metabaser)
metabase.connect(dbname = "prd380", uid = "mb_user", pwd = "Temp_usr_9999", host = "metabase.am.lilly.com", port = 1525, type = "sid", driver = "jdbc")


upstream_interactions <-
function(gene, distance)
{
	#create an empty frame that will be holding the result
	result <- data.frame(
	gene1 = character(),
	gene2 = character(),
	distance = numeric(),
	direction = character(),
	effect = character(),
	mechanisim = character(),
	trust = character()
	)
	gene_list <- c(gene)
	counter <- 0
	while (counter < distance)
	{
		counter <- counter + 1
		gene_list_temp <- c(character(0))
		for (gene in gene_list)
		{
			if(length(recognize.gene.symbols(gene, drop.unrecognized = TRUE)) > 0)
			{
				gene_id <- recognize.gene.symbols(toupper(gene))
				upstream <- get.gene.interactions(gene_id, filter = list(effect = c("Activation", "Inhibition"), has.compounds = FALSE, has.genes = TRUE), direction = 'upstream')
				fr <- data.frame(
				gene1 = upstream$id1,
				gene2 = upstream$id2,
				distance = rep(counter, length(upstream$id1)),
				direction = rep("upstream", length(upstream$id1)),
				effect = upstream$effect,
				mechanisim = upstream$mechanism,
				trust = upstream$trust
				)
				result <- rbind(result, fr)
				gene_list_temp <- c(gene_list_temp, upstream$id1)
			}
		}
		gene_list <- gene_list_temp
	}
	return(result)
}

show_interactions <-
function(gene, distance, direct)
{
	#create an empty frame that will be holding the result
	result <- data.frame(
	gene1 = character(),
	gene2 = character(),
	distance = numeric(),
	direction = character(),
	effect = character(),
	mechanisim = character(),
	trust = character()
	)
	gene_list <- c(gene)
	counter <- 0
	while (counter < distance)
	{
		counter <- counter + 1
		gene_list_temp <- c(character(0))
		for (gene in gene_list)
		{
			if(length(recognize.gene.symbols(toupper(gene))) > 0)
			{
				gene_id <- recognize.gene.symbols(gene)
				stream <- get.gene.interactions(gene_id, filter = list(effect = c("Activation", "Inhibition"), has.compounds = FALSE, has.genes = TRUE), direction = direct)
				fr <- data.frame(
				gene1 = stream$id1,
				gene2 = stream$id2,
				distance = rep(counter, length(stream$id1)),
				direction = rep(direct, length(stream$id1)),
				effect = stream$effect,
				mechanisim = stream$mechanism,
				trust = stream$trust
				)
				result <- rbind(result, fr)
				if(toupper(direct) == "UPSTREAM")
					gene_list_temp <- c(gene_list_temp, stream$id1)
				else
					gene_list_temp <- c(gene_list_temp, stream$id2)
				
			}
		}
		gene_list <- gene_list_temp
	}
	return(result)
}
