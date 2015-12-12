GENERAL-INFO-START

	seq-file            seqs-sample.txt
	trace-file          logs\sample-data.trace.tsv		
	coal-stats-file		logs\sample-data.flatStats.tsv
#	num-pop-partitions	1
	locus-mut-rate		CONST

	mcmc-iterations	  200
	iterations-per-log  5
	logs-per-line       10


	find-finetunes		FALSE
	finetune-coal-time	0.01		
	finetune-mig-time	0.3		
	finetune-theta		0.04
	finetune-mig-rate	0.02
	finetune-tau		0.0000008
	finetune-mixing		0.003
#   finetune-locus-rate 0.3
	
	tau-theta-print		10000.0
	tau-theta-alpha		1.0			# for STD/mean ratio of 100%
	tau-theta-beta		10000.0		# for mean of 1e-4

	mig-rate-print		0.001
	mig-rate-alpha		0.002
	mig-rate-beta		0.00001

GENERAL-INFO-END

CURRENT-POPS-START	

	POP-START
		name		A
		samples		one d
	POP-END

	POP-START
		name		B
		samples		two d
	POP-END

	POP-START
		name		C
		samples		three d
	POP-END

	POP-START
		name		D
		samples		five d
	POP-END
	
CURRENT-POPS-END

#!! ROOT population must be placed last!
ANCESTRAL-POPS-START

	POP-START
		name			AB
		children		A		B
		tau-initial	0.000005
		tau-beta		20000.0	
		finetune-tau			0.0000008
	POP-END

	POP-START
		name			ABC
		children		AB		C
		tau-initial	0.00001
		tau-beta		20000.0	
		finetune-tau			0.0000008
	POP-END

	POP-START
		name			root
		children		ABC	D
		tau-initial	0.00005
		tau-beta		20000.0	
		finetune-tau			0.00000286
	POP-END

ANCESTRAL-POPS-END


MIG-BANDS-START	

#	BAND-START		
#       source  D
#       target  B
#       mig-rate-print 0.1
#	BAND-END

MIG-BANDS-END
