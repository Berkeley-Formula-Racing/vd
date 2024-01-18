function torque_fn = KTM450()

% Stock
points = [
   0 0           
% 19.95907413	4632.81866
% 19.95192271	4632.471246
% 19.94700219	4632.420599
% 19.94445671	4632.677193
% 19.94420227	4633.236416
% 19.94679543	4634.138664
% 19.95212314	4635.359485
% 19.96038728	4636.911536
% 19.97097499	4638.728066
% 19.98395398	4640.784973
% 19.99891981	4643.053187
% 20.01583665	4645.522985
% 20.03436536	4648.154108
% 20.05450332	4650.931921
% 20.07564785	4653.785806
% 20.09766318	4656.710493
% 20.12029993	4659.680247
% 20.14328117	4662.667765
% 20.16647744	4665.661576
% 20.19008446	4668.694255
% 20.21359503	4671.702597
% 20.23685386	4674.671186
% 20.26018551	4677.643553
% 20.282994	4680.551195
% 20.30510155	4683.380618
% 20.32631532	4686.118953
% 20.34644888	4688.752201
% 20.36495696	4691.217581
% 20.38251545	4693.595183
% 20.39897235	4695.865707
% 20.41447039	4698.035539
% 20.42924505	4700.133162
% 20.44341403	4702.174538
% 20.45681816	4704.140224
% 20.469729	4706.066308
% 20.48192961	4707.926338
% 20.4938193	4709.774476
% 20.50579638	4711.671446
% 20.51776159	4713.605148
% 20.52984329	4715.596596
% 20.54247376	4717.71725
% 20.55536032	4719.927267
% 20.56817612	4722.17595
% 20.58129448	4724.518542
% 20.59473427	4726.955675
% 20.60853425	4729.485348
% 20.62269582	4732.098482
% 20.6372544	4734.792109
% 20.65222578	4737.562292
% 20.66759168	4740.403061
% 20.68360712	4743.351173
% 20.69947635	4746.307867
% 20.71638129	4749.40108
% 20.73275502	4752.41467
% 20.74883054	4755.391957
% 20.76441553	4758.289009
% 20.7799867	4761.162115
% 20.79465816	4763.922947
% 20.80926388	4766.6926
% 20.82329889	4769.3975
% 20.83647499	4772.009937
% 20.8487203	4774.521111
% 20.85991067	4776.904855
% 20.87015978	4779.180111
% 20.87952868	4781.349236
% 20.88794269	4783.401307
% 20.89548659	4785.346163
% 20.90267181	4787.240637
% 20.9088488	4789.001267
% 20.91427234	4790.65729
% 20.91929218	4792.246001
% 20.92386395	4793.762965
% 20.92752882	4795.167031
% 20.93129355	4796.56545
% 20.93478597	4797.941838
% 20.93811068	4799.291472
% 20.94127331	4800.633583
% 20.94438618	4801.987297
% 20.94758484	4803.375429
% 20.95083735	4804.781378
% 20.95423411	4806.241146
% 20.95789827	4807.755393
% 20.96197107	4809.344281
% 20.96659689	4811.015817
% 20.97188497	4812.769486
% 20.97791789	4814.613146
% 20.98470593	4816.546076
% 20.99218581	4818.555214
% 21.00022488	4820.626261
% 21.00889019	4822.778385
% 21.01808513	4825.003982
% 21.02772556	4827.297237
% 21.03767785	4829.646424
% 21.04793527	4832.053864
% 21.05830355	4834.499167
% 21.06870366	4836.972026
% 21.07964861	4839.584962
% 21.08983508	4842.071656
% 21.0997935	4844.543822
% 21.10973644	4847.050924
% 21.11917634	4849.507407
% 21.12751473	4851.78476
% 21.1358213	4854.116045
% 21.14340634	4856.364926
% 21.14992398	4858.463712
% 21.15569773	4860.47287
% 21.16053409	4862.374091
% 21.16438287	4864.163209
% 21.16717523	4865.8368
% 21.1690165	4867.434795
% 21.16952655	4868.870977
% 21.16888616	4870.194622
% 21.16732558	4871.52926
% 21.16421047	4872.657294
% 21.15987487	4873.667873
% 21.15439179	4874.623482
% 21.14769293	4875.506571
% 21.13979315	4876.263662
% 21.13094348	4877.021606
% 21.12092937	4877.730167
% 21.10979183	4878.395924
% 21.09773117	4879.020654
% 21.08501033	4879.606919
% 21.07105836	4880.189272
% 21.05646218	4880.754399
% 21.04136853	4881.308301
% 21.02578489	4881.862516
% 21.00996625	4882.418975
% 20.99430501	4882.96709
% 20.97855333	4883.524165
% 20.96272248	4884.094801
% 20.9468003	4884.682653
% 20.93087225	4885.26719
% 20.91456292	4885.89678
% 20.89826336	4886.546254
% 20.88202789	4887.233378
% 20.86582114	4887.950176
% 20.849699	4888.722164
% 20.83410683	4889.510033
% 20.81877075	4890.358544
% 20.80377176	4891.260888
% 20.78926307	4892.233379
% 20.77520549	4893.286934
% 20.76164409	4894.416268
% 20.74866936	4895.626235
% 20.7360767	4896.947701
% 20.72423688	4898.364242
% 20.71325374	4899.886099
% 20.70323903	4901.525843
% 20.69433814	4903.295718
% 20.68702393	4905.187815
% 20.68106337	4907.237526
% 20.67664414	4909.451707
% 20.67385372	4911.838473
% 20.67273116	4914.405393
% 20.67329343	4917.159777
% 20.67554494	4920.097639
% 20.67969933	4923.277414
% 20.6853772	4926.667422
% 20.69262357	4930.260867
% 20.70147041	4934.044383
% 20.71195617	4938.020321
% 20.72322578	4942.318937
% 20.73652402	4946.578865
% 20.75111145	4950.978461
% 20.76684068	4955.520239
% 20.78355972	4960.190514
% 20.80072873	4964.68665
% 20.81899212	4969.572991
% 20.8379804	4974.565955
% 20.85759807	4979.657403
% 20.87775814	4984.840632
% 20.89840231	4990.110776
% 20.91949834	4995.465055
% 20.94105112	5000.897849
% 20.96322026	5006.420621
% 20.98591625	5012.013425
% 21.00962122	5017.786355
% 21.03400168	5023.628926
% 21.05909566	5029.543307
% 21.08455028	5035.445403
% 21.11092828	5041.471171
% 21.1379158	5047.542649
% 21.16549051	5053.6665
% 21.19360958	5059.840469
% 21.22259879	5066.147944
% 21.25147802	5072.402626
% 21.28002406	5078.559625
% 21.30850675	5084.698261
% 21.33691532	5090.813831
% 21.36534163	5096.911918
% 21.39355523	5102.957354
% 21.42150272	5108.937135
% 21.44903178	5114.826615
% 21.47604169	5120.607858
% 21.50240264	5126.262765
% 21.52795705	5131.772037
% 21.55223195	5137.077298
% 21.57531568	5142.198861
% 21.59704973	5147.123953
% 21.61732407	5151.846664
% 21.63596275	5156.34774
% 21.65389013	5160.804062
% 21.66987533	5165.009645
% 21.68382966	5168.959971
% 21.69592754	5172.685354
% 21.70566669	5176.137381
% 21.71281971	5179.25086
% 21.71786172	5182.143105
% 21.72078011	5184.824901
% 21.72138522	5187.283462
% 21.72014846	5189.577954
% 21.71699465	5191.699662
% 21.71214496	5193.678603
% 21.70578094	5195.52989
% 21.69800665	5197.257117
% 21.68927499	5198.904433
% 21.67971259	5200.478813
% 21.66953005	5201.991331
% 21.65892259	5203.459297
% 21.64786669	5204.931491
% 21.63658065	5206.374208
% 21.62516742	5207.798696
% 21.61363613	5209.213737
% 21.6019759	5210.628588
% 21.59031463	5212.030839
% 21.57882129	5213.406612
% 21.56693434	5214.814099
% 21.55529182	5216.201787
% 21.54364584	5217.578756
% 21.53188925	5218.975599
% 21.51994779	5220.404031
% 21.50834827	5221.818082
% 21.49663321	5223.264884
% 21.48504716	5224.73779
% 21.47388322	5226.201891
% 21.46286283	5227.689344
% 21.4520175	5229.188695
% 21.44123819	5230.722853
% 21.43051544	5232.274501
% 21.41977491	5233.830696
% 21.4089703	5235.400148
% 21.39782275	5237.01437
% 21.38649735	5238.628089
% 21.37505196	5240.259313
% 21.36356633	5241.922642
% 21.35199434	5243.609904
% 21.34066844	5245.298096
% 21.32921229	5247.021758
% 21.31767196	5248.785753
% 21.30648001	5250.573719
% 21.2950916	5252.435948
% 21.28394148	5254.343923
% 21.27320379	5256.293171
% 21.26294105	5258.282965
% 21.25293063	5260.352567
% 21.24373194	5262.42729
% 21.23511204	5264.54821
% 21.22709998	5266.710052
% 21.21969692	5268.916589
% 21.21269668	5271.222219
% 21.20636768	5273.534361
% 21.20048373	5275.888993
% 21.19491455	5278.291162
% 21.1895732	5280.730297
% 21.18455916	5283.133333
% 21.17956285	5285.59039
% 21.17469074	5288.049058
% 21.16990656	5290.502736
% 21.16503188	5293.003392
% 21.16013059	5295.490918
% 21.15514414	5297.960429
% 21.15002647	5300.406106
% 21.14465326	5302.87248
% 21.13928996	5305.187601
% 21.13360649	5307.489618
% 21.12758578	5309.75674
% 21.12112792	5311.990106
% 21.11429577	5314.144203
% 21.10663811	5316.321454
% 21.09867646	5318.451608
% 21.09019265	5320.552512
% 21.08132567	5322.601087
% 21.07177971	5324.695746
% 21.06203174	5326.771404
% 21.0517635	5328.789881
% 21.04106925	5330.800153
% 21.02991873	5332.829869
% 21.01878632	5334.792293
% 21.00737521	5336.75754
% 20.99557177	5338.77176
% 20.9836273	5340.801183
% 20.97165376	5342.853828
% 20.9597088	5344.936313
% 20.94786108	5347.051619
% 20.93587883	5349.246522
% 20.92399554	5351.473022
% 20.91212358	5353.731327
% 20.90014128	5356.028408
% 20.88800014	5358.351855
% 20.87581873	5360.651525
% 20.8632831	5362.973178
% 20.85029627	5365.320259
% 20.83687062	5367.681011
% 20.82287245	5370.07906
% 20.80834885	5372.509962
% 20.79312624	5375.008401
% 20.77749919	5377.537726
% 20.76182372	5380.0573
% 20.7463487	5382.556752
% 20.73086403	5385.102768
% 20.71572924	5387.669489
% 20.70053032	5390.363024
% 20.68544232	5393.197379
% 20.67066379	5396.186292
% 20.65674103	5399.281887
% 20.64385954	5402.490471
% 20.63251274	5405.752398
% 20.62258139	5409.140756
20.61404568	5412.646093
20.60690442	5416.264469
20.60111742	5419.989656
20.59644008	5423.860724
20.59279767	5427.805443
20.58947433	5432.055532
20.58799849	5436.174306
20.58755421	5440.376796
20.58827279	5444.612513
20.58997553	5448.9352
20.59335952	5453.06434
20.59740082	5457.541986
20.6023849	5462.162703
20.60845872	5466.893939
20.61574877	5471.666699
20.62402604	5476.61541
20.63332536	5481.698566
20.64421697	5487.036216
20.65657747	5492.594194
20.67059132	5498.501376
20.68643023	5504.709651
20.70406659	5511.231237
20.72315733	5517.940404
20.74415693	5524.995601
20.76716913	5532.417693
20.79229203	5540.224191
20.81960517	5548.430461
20.84919394	5557.05124
20.8811496	5566.097127
20.9157915	5575.613953
20.95333229	5585.624577
20.99331023	5596.031139
21.03677661	5607.075053
21.08330055	5618.576421
21.13284914	5630.524062
21.18547337	5642.913306
21.24211487	5655.885141
21.30159979	5669.150029
21.36465438	5682.909922
21.43168927	5697.209791
21.50306596	5712.101106
21.57772383	5727.416923
21.65616881	5743.241188
21.73861943	5759.614381
21.8259245	5776.750801
21.91646461	5794.325055
22.01126917	5812.509131
22.10970929	5831.211907
22.21170361	5850.430769
22.31583312	5869.893062
22.42310272	5889.834367
22.53306887	5910.187766
22.64589333	5931.00157
22.76138825	5952.258642
22.8788812	5973.863513
22.99886696	5995.912863
23.12395072	6018.912046
23.2512059	6042.334002
23.38064767	6066.191008
23.51223144	6090.481321
23.64433453	6114.906929
23.77553526	6139.200827
23.91128443	6164.387657
24.0487145	6189.931807
24.18775202	6215.81694
24.33001638	6242.342602
24.4738599	6269.197597
24.61614357	6295.785312
24.75972719	6322.644868
24.90450282	6349.75501
25.05027665	6377.080993
25.19702243	6404.608606
25.34460374	6432.309394
25.49286572	6460.152494
25.64379669	6488.506832
25.79773474	6517.425572
25.94854242	6545.763032
26.09931529	6574.097252
26.24963288	6602.369085
26.39689491	6630.114023
26.54035257	6657.211341
26.68606648	6684.773902
26.83014963	6712.110787
26.97078825	6738.887046
27.11092921	6765.658957
27.25229362	6792.758673
27.38741809	6818.856042
27.5195868	6844.564523
27.65032298	6870.17676
27.7765114	6895.118429
27.89620814	6919.032401
28.01607547	6943.15185
28.13284549	6966.871377
28.24652608	6990.194607
28.35988858	7013.673336
28.46978493	7036.717723
28.5762221	7059.331707
28.67920881	7081.519633
28.77855596	7103.247356
28.87150271	7123.949431
28.96155502	7144.329487
29.04926842	7164.510219
29.13320002	7184.195228
29.21398631	7203.515935
29.2913031	7222.448727
29.36484511	7240.875642
29.43428502	7258.708245
29.50194374	7276.47051
29.5663049	7293.811522
29.62785054	7310.791233
29.68693876	7327.541514
29.74285675	7343.880076
29.79483939	7359.588512
29.84390864	7374.907685
29.890138	7389.835201
29.93103022	7403.923424
29.97164472	7417.988506
30.00949803	7431.652733
30.0453453	7445.027071
30.07890717	7458.07072
30.10916721	7470.632568
30.13790362	7482.981281
30.16530853	7495.149156
30.19096204	7507.192828
30.21480792	7519.018238
30.2385461	7530.845772
30.25904502	7542.172136
30.27774073	7553.205243
30.29467808	7563.821961
30.3105324	7574.204894
30.32526724	7584.384964
30.33921339	7594.435252
30.35228521	7604.386906
30.3647409	7614.29939
30.37624681	7624.068057
30.38695734	7633.725102
30.3967544	7643.252185
30.40571176	7652.6297
30.41383472	7661.855988
30.42127452	7671.02363
30.42822833	7680.149121
30.43448231	7689.230131
30.44027462	7698.278647
30.4454236	7707.300614
30.44988888	7716.284995
30.45346627	7725.236011
30.45620198	7734.156102
30.45793222	7743.042429
30.4587158	7751.88494
30.45866327	7760.901612
30.45764567	7769.882602
30.45580086	7778.831268
30.45294383	7787.734759
30.44928841	7796.604007
30.44471332	7805.226422
30.43929513	7813.799003
30.43307105	7822.318171
30.42605586	7830.779034
30.41822607	7839.174902
30.40976179	7847.514909
30.40016434	7856.953716
30.39044838	7865.171912
30.38026585	7873.53785
30.36962761	7881.652252
30.35867173	7889.726285
30.34667951	7898.49849
30.33462668	7906.853393
30.32272474	7914.611198
30.31043349	7922.525637
30.29775564	7930.400196
30.28689429	7936.439193
30.27416384	7943.911306
30.2606434	7951.720068
30.24666293	7959.699901
30.23259252	7967.484366
30.21803924	7975.441401
30.20356584	7983.220546
30.18870943	7991.181687
30.17412259	7998.893221
30.15924452	8006.772867
30.14449957	8014.582763
30.12943797	8022.590921
30.11444949	8030.604356
30.09941513	8038.70416
30.08437647	8046.878802
30.0702331	8054.761873
30.05602808	8062.821337
30.04287266	8070.560757
30.03009983	8078.320665
30.01774183	8086.05349
30.00523232	8094.019389
29.99315926	8102.010387
29.98119797	8110.203568
29.96973355	8118.429347
29.95892919	8126.698677
29.94869951	8135.006711
29.93941904	8143.181267
29.93053035	8151.518619
29.92278624	8159.730727
29.91558014	8167.971831
29.90907978	8176.251135
29.90329738	8184.566872
29.89815968	8192.967894
29.89335675	8201.569288
29.88944761	8210.212476
29.88627572	8218.884286
29.88375355	8227.669742
29.88218167	8236.243728
29.88115972	8245.028419
29.88093716	8253.659416
29.8814106	8262.327428
29.88238064	8270.933176
29.88378105	8279.629677
29.88558042	8288.162621
29.88766275	8296.910177
29.88998415	8305.677522
29.89262964	8314.464617
29.89556713	8323.456621
29.89874412	8332.255438
29.90224892	8341.244086
29.90604145	8350.214452
29.9101051	8359.158484
29.91431059	8367.862346
29.9188006	8376.717762
29.92332996	8385.364636
29.92789282	8393.908283
29.93246755	8402.388741
29.93702413	8410.805328
29.94154621	8419.154814
29.9459924	8427.388586
29.95053438	8435.781247
29.95499729	8444.095474
29.9593471	8452.302977
29.96345997	8460.406636
29.96717436	8468.231681
29.97038513	8475.77446
29.97312877	8483.21769
29.97510517	8490.540698
29.97671002	8497.813097
29.97772275	8505.152621
29.97797068	8512.382143
29.97740695	8519.50058
29.97624482	8526.549158
29.97401305	8533.624527
29.97061923	8540.400503
29.96641435	8547.093242
29.96126585	8553.685038
29.95518513	8560.176274
29.94808936	8566.39224
29.94025191	8572.689309
29.93144914	8578.873124
29.92187065	8584.960316
29.9115833	8590.951865
29.9007718	8596.863035
29.88939967	8602.678532
29.87737517	8608.542324
29.86503528	8614.239085
29.85236524	8619.864271
29.83945813	8625.390539
29.82646101	8630.822991
29.81371352	8636.031264
29.80077767	8641.273632
29.78805073	8646.414813
29.7755613	8651.482339
29.7633282	8656.486945
29.75134925	8661.42456
29.73977953	8666.245075
29.72847792	8671.011302
29.71746749	8675.722205
29.70682231	8680.380316
29.69655253	8684.999235
29.68674896	8689.582674
29.67746278	8694.139052
29.66865623	8698.679683
29.66034303	8703.210494
29.65253731	8707.737654
29.64515322	8712.26612
29.63850028	8716.724687
29.63206705	8721.313589
29.62623028	8725.886261
29.62088802	8730.464005
29.61606055	8735.057394
29.61074515	8739.930672
29.60630589	8744.502287
29.60188714	8749.120712
29.59734498	8753.77051
29.59244239	8758.444496
29.5876899	8762.962107
29.58207842	8767.720001
29.57572069	8772.624664
29.56903245	8777.483091
29.56181718	8782.403532
29.55407728	8787.391352
29.54563789	8792.555547
29.53691842	8797.690778
29.52760347	8803.009522
29.51797118	8808.406548
29.50812711	8813.87974
29.49835432	8819.315028
29.48859373	8824.816166
29.47887549	8830.39131
29.46934776	8836.024963
29.46008459	8841.714595
29.45108075	8847.523264
29.44244933	8853.39133
29.43432012	8859.304052
29.42656166	8865.387279
29.41927035	8871.521793
29.41249546	8877.633553
29.40621515	8883.652064
29.40034023	8889.690983
29.39483927	8895.61921
29.38938508	8901.441548
29.38428248	8907.233394
29.37914245	8913.139425
29.37385699	8919.023421
29.36825087	8925.007486
29.36238562	8931.087638
29.35580906	8937.133339
29.34865224	8943.134924
29.3408612	8949.088315
29.33253786	8954.85684
29.32354661	8960.575014
29.31388721	8966.229045
29.30351939	8971.809804
29.29252771	8977.320602
29.28072342	8982.894068
29.26863297	8988.249301
29.25609141	8993.528464
29.24294619	8998.838409
29.22931647	9004.103152
29.21567791	9009.142447
29.2014412	9014.210358
29.18697622	9019.178039
29.17263234	9023.951617
29.15835637	9028.58737
29.14418757	9033.109428
29.13003407	9037.571849
29.11607637	9041.949469
29.10239796	9046.245555
29.08907686	9050.463064
29.07609185	9054.623062
29.06373781	9058.683827
29.0519982	9062.666826
29.04098022	9066.570151
29.03046388	9070.493283
29.02104703	9074.306007
29.01271745	9078.023653
29.005555	9081.64273
28.99961632	9085.159098
28.99512308	9088.482058
28.99171643	9091.72466
28.989481	9094.884286
28.98843513	9097.945025
28.98853929	9100.920406
28.98977983	9103.78735
28.99209529	9106.541824
28.99547009	9109.373873
28.9997503	9111.888369
29.00480684	9114.287049
29.01053065	9116.589793
29.01684076	9118.793443
29.02325702	9120.754181
29.03015058	9122.788859
29.03712541	9124.73961
29.04406701	9126.611742
29.05081781	9128.413571
29.05740458	9130.149415
29.06382196	9131.827123
29.07014077	9133.491151
29.07633405	9135.110197
29.0824341	9136.690306
29.0884779	9138.237957
29.09450125	9139.757785
29.10048887	9141.218872
29.10659577	9142.655516
29.11278336	9144.072798
29.11911742	9145.484586
29.12565213	9146.882007
29.13240095	9148.270538
29.1393804	9149.659132
29.14680272	9151.079467
29.15433552	9152.455228
29.16210794	9153.838605
29.17005679	9155.232088
29.17808627	9156.639986
29.18589288	9158.034044
29.19382851	9159.487814
29.20161555	9160.96625
29.20920229	9162.468022
29.21654801	9163.989649
29.2238698	9165.563971
29.23081338	9167.147228
29.23702278	9168.727768
29.24293618	9170.31148
29.24865464	9171.89398
29.25374352	9173.432162
29.25850253	9174.959712
29.26321897	9176.474273
29.26744053	9177.972481
29.27108253	9179.455162
29.2743249	9180.920464
29.27712873	9182.364813
29.27944403	9183.783635
29.28121349	9185.17151
29.28239101	9186.522401
29.28303923	9187.858203
29.2829909	9189.144102
29.28221303	9190.375658
29.28069153	9191.54999
29.27839952	9192.667097
29.27516374	9193.704409
29.27116732	9194.692779
29.26640726	9195.63786
29.260885	9196.546514
29.25461255	9197.426183
29.24770491	9198.286199
29.24013055	9199.131249
29.23196401	9199.968294
29.22324698	9200.800268
29.21401362	9201.658351
29.20443557	9202.522167
29.19463614	9203.392643
29.18460273	9204.283427
29.17452704	9205.188184
29.1647571	9206.076621
29.1551735	9206.972943
29.14584713	9207.876132
29.13698568	9208.770724
29.12852623	9209.66474
29.12048936	9210.556571
29.11289353	9211.444789
29.10572825	9212.328123
29.09897585	9213.204838
29.09262332	9214.073562
29.08665365	9214.932113
29.08104822	9215.777781
29.07567013	9216.62391
29.0700323	9217.537243
29.06492119	9218.400877
29.06041708	9219.185897
29.05624674	9219.924157
29.05233407	9220.614409
29.04908154	9221.191311
29.04587727	9221.751445
29.04260067	9222.30738
29.03939711	9222.825895
29.03640363	9223.293718
29.03355633	9223.716587
29.03083416	9224.096029
29.02822481	9224.433192
29.0256983	9224.731382
29.02322165	9224.99422
29.02074997	9225.226293
29.01825863	9225.432111
29.01573304	9225.617351
29.01309434	9225.789548
29.01044699	9225.95153
29.00782363	9226.109607
29.00524929	9226.270011
29.00273559	9226.438048
29.00033604	9226.617902
28.99797808	9226.813923
28.99562813	9227.029474
28.99319697	9227.2722
28.99071082	9227.538031
28.98806879	9227.83309
28.98522862	9228.159623
28.98214787	9228.520429
28.97885461	9228.907802
28.97522884	9229.339274
28.9712884	9229.80632
28.9670067	9230.311057
28.96239335	9230.852669
28.95743609	9231.433421
28.95211371	9232.054609
28.94641244	9232.718857
28.94028714	9233.428776
28.93366037	9234.189409
28.92634307	9235.010348
28.91844734	9235.89167
28.90991739	9236.83873
28.9007074	9237.856482
28.89078125	9238.948475
28.88022851	9240.108345
28.86891229	9241.347143
28.85657102	9242.697513
28.84363791	9244.116514
28.83005274	9245.622584
28.81583134	9247.231851
28.80122648	9248.913582
28.78660594	9250.66399
28.77143894	9252.520399
28.75599748	9254.461827
28.74033972	9256.471811
28.72441159	9258.582478
28.70830571	9260.757461
28.6919896	9263.012203
28.67547534	9265.344767
28.65879272	9267.759869
28.64199508	9270.241138
28.62515742	9272.770883
28.60833728	9275.360408
28.59164052	9278.001546
28.57509619	9280.6845
28.55873044	9283.417396
28.54257978	9286.216666
28.5265849	9289.056611
28.51065837	9291.941504
28.49450575	9294.906876
28.47830079	9297.924134
28.46189213	9301.00152
28.44542635	9304.158417
28.42871161	9307.369915
28.41211597	9310.608167
28.39535995	9313.917817
28.3785534	9317.302514
28.36154287	9320.747303
28.34461506	9324.26043
28.32763452	9327.831927
28.31051585	9331.468375
28.29327545	9335.147692
28.27589411	9338.875386
28.25827654	9342.670448
28.2400002	9346.612572
28.22155798	9350.588304
28.20288855	9354.601141
28.18391434	9358.643475
28.16466029	9362.699883
28.14554141	9366.6753
28.12592979	9370.638852
28.10585502	9374.581641
28.08528214	9378.489636
28.06414731	9382.353423
28.04192407	9386.230078
28.01959911	9389.946735
27.99642613	9393.681534
27.97291726	9397.338857
27.9491512	9400.917064
27.92569854	9404.346659
27.90169664	9407.785483
27.87800003	9411.06139
27.85418294	9414.26327
27.83017637	9417.400726
27.80562679	9420.531873
27.78082668	9423.59603
27.75609327	9426.529314
27.73103841	9429.413569
27.70567038	9432.250586
27.68037425	9435.006212
27.65473698	9437.750834
27.62781544	9440.618117
27.60047315	9443.487638
27.5726954	9446.369239
27.54443903	9449.269148
27.51572138	9452.191624
27.48710542	9455.082756
27.45688066	9458.134852
27.42749354	9461.094166
27.39783702	9464.085794
27.36803574	9467.10088
27.3373022	9470.222071
27.30761611	9473.253781
27.27635812	9476.46394
27.2448082	9479.73127
27.21305412	9483.062947
27.18204572	9486.374313
27.15148274	9489.703281
27.1206112	9493.128255
27.08898701	9496.693202
27.05865461	9500.200331
27.0283062	9503.763491
26.99771678	9507.415778
26.96809658	9511.046878
26.94070899	9514.538072
26.91261951	9518.155109
26.88529327	9521.779788
26.85871293	9525.425271
26.83262732	9529.121692
26.805681	9533.024116
26.77939656	9536.998687
26.75404041	9541.094938
26.72936265	9545.263367
26.7055888	9549.55288
26.68323177	9553.872051
26.66162283	9558.299828
26.64061729	9562.821855
26.62035654	9567.487326
26.60059511	9572.261283
26.58137315	9577.173769
26.56263718	9582.233477
26.54408814	9587.496576
26.52583184	9592.947138
26.50774984	9598.614912
26.4900478	9604.456054
26.47276292	9610.481802
26.45608669	9616.615513
26.43991168	9622.874122
26.42432404	9629.21888
26.40910994	9635.686632
26.39424277	9642.26701
26.37975111	9648.981916
26.36559308	9655.810994
26.35178103	9662.749041
26.33838857	9669.805888
26.32543861	9676.970349
26.31287901	9684.22177
26.30073601	9691.570917
26.28897981	9699.012382
26.27755272	9706.533559
26.26613002	9714.282421
26.2550223	9722.092095
26.24422992	9729.946803
26.23373982	9737.83092
26.22353285	9745.728907
26.21362477	9753.633718
26.20414144	9761.386003
26.19485804	9769.102222
26.18575046	9776.786724
26.1768699	9784.381999
26.1683661	9791.719066
26.16004466	9798.914763
26.15180424	9805.986372
26.14361371	9812.898894
26.13538024	9819.682507
26.12707056	9826.315125
26.11830043	9833.142221
26.10949121	9839.62359
26.10023072	9846.111285
26.09060495	9852.43027
26.08052906	9858.574635
26.07013906	9864.398267
26.05901943	9870.255996
26.04748671	9875.804967
26.03530593	9881.218997
26.02244704	9886.498706
26.00890906	9891.64518
25.99473645	9896.661559
25.98002392	9901.554737
25.96487844	9906.331824
25.94913906	9911.058229
25.93316167	9915.723316
25.91693569	9920.289207
25.90046353	9924.762138
25.88383631	9929.06746
25.86728836	9933.24475
25.85012843	9937.415391
25.83292306	9941.542166
25.81565576	9945.591785
25.79816821	9949.66374
25.78045521	9953.678852
25.76298058	9957.53898
25.7452727	9961.329244
25.72736045	9965.072297
25.70933282	9968.7391
25.69123421	9972.333175
25.67304682	9975.855104
25.65478194	9979.303656
25.6364646	9982.682688
25.61809363	9985.995647
25.59927062	9989.313749
25.58027868	9992.570993
25.56107255	9995.768304
25.54160666	9998.909383
25.52185489	10001.99833
25.5023673	10004.96625
25.48227494	10007.94169
25.46180374	10010.88281
25.44088731	10013.79249
25.41944459	10016.67294
25.39737908	10019.52629
25.37497077	10022.30787
25.35146501	10025.10201
25.32727613	10027.86774
25.30222829	10030.62172
25.27610563	10033.39723
25.2491613	10036.14785
25.22182018	10038.8469
25.19356687	10041.54509
25.16461885	10044.23383
25.13513545	10046.9201
25.10501052	10049.63121
25.07439266	10052.37551
25.0433878	10055.16124
25.01203673	10057.99806
24.97978001	10060.96525
24.94788801	10063.92977
24.9158079	10066.96347
24.88359923	10070.06832
24.8507472	10073.29512
24.8186617	10076.50339
24.78526797	10079.93798
24.7524961	10083.32021
24.71990493	10086.74577
24.688077	10090.12189
24.65610367	10093.55776
24.6250447	10096.93489
24.59391237	10100.40103
24.56335792	10103.83683
24.53345827	10107.24805
24.50459456	10110.57637
24.47661483	10113.84001
24.44929704	10117.06295
24.42262167	10120.24129
24.39649608	10123.37573
24.370867	10126.46183
24.34563559	10129.50022
24.32072108	10132.49188
24.29606163	10135.43862
24.27162575	10138.34218
24.24741241	10141.20268
24.22279348	10144.08131
24.19837108	10146.91561
24.17412285	10149.70639
24.15001579	10152.45727
24.12610903	10155.16757
24.10307243	10157.79146
24.08042587	10160.39804
24.05808269	10162.99012
24.03611074	10165.57629
24.01354629	10168.24693
23.99225301	10170.89071
23.97133797	10173.54821
23.95110608	10176.23537
23.93166198	10178.95508
23.91390086	10181.64526
23.89640769	10184.41574
23.88005323	10187.22863
23.86474119	10190.07868
23.85068044	10192.97502
23.83807531	10195.91017
23.82673128	10198.89151
23.81666362	10201.9197
23.80796297	10204.99732
23.80019219	10208.18237
23.79393604	10211.34043
23.7889421	10214.53399
23.78498562	10217.83029
23.78242278	10221.09587
23.7809764	10224.38146
23.78045544	10227.75961
23.78090517	10231.14812
23.78223135	10234.45463
23.78424007	10237.82089
23.78683761	10241.12505
23.78988536	10244.38256
23.79325508	10247.57686
23.79683028	10250.69264
23.80049825	10253.71009
23.80409671	10256.58252
23.80758398	10259.35422
23.81086681	10262.0226
23.81386074	10264.58799
23.81648997	10267.05245
23.81869605	10269.41953
23.82041691	10271.69462
23.82156033	10273.87351
23.8221835	10275.98652
23.82195725	10277.97567
23.82141826	10279.96079
23.82025649	10281.88083
23.81858999	10283.75693
23.81635204	10285.62493
23.81397235	10287.48951
23.81079045	10289.27052
23.80720847	10291.01794
23.80329176	10292.77767
23.79899442	10294.44868
23.79436005	10296.09461
23.78921791	10297.70069
23.78375661	10299.29895
23.77772846	10300.82924
23.77127691	10302.35755
23.7644065	10303.84391
23.7572886	10305.33051
23.74978071	10306.81823
23.7421642	10308.29606
23.73440462	10309.762
23.72648809	10311.24264
23.71848734	10312.71843
23.71041868	10314.17512
23.70222371	10315.63654
23.69389343	10317.10334
23.68544242	10318.57827
23.67687675	10320.05897
23.66821642	10321.54087
23.65943237	10323.02269
23.65048151	10324.5056
23.64133632	10325.99093
23.63196872	10327.48309
23.62233953	10328.9857
23.61239334	10330.5023
23.60208688	10332.03372
23.59136063	10333.58073
23.58015514	10335.14305
23.56842423	10336.72039
23.55612993	10338.31343
23.54323474	10339.92377
23.52973328	10341.55374
23.51565068	10343.20684
23.50103312	10344.8873
23.48561711	10346.63195
23.4701651	10348.37625
23.45234251	10350.36842
23.43637006	10352.18196
23.41991568	10354.05495
23.40373077	10355.90664
23.38715217	10357.80764
23.37265922	10359.4718
23.35609308	10361.37691
23.33968749	10363.26724
23.32370429	10365.1115
23.3078481	10366.95078
23.29249724	10368.74287
23.27702817	10370.56976
23.26243753	10372.32005
23.24752704	10374.15704
23.23314013	10375.98993
23.21908021	10377.862
23.2060803	10379.68949
23.19358641	10381.55868
23.18225188	10383.38123
23.1717783	10385.20265
23.16207638	10387.0285
23.15322804	10388.84949
23.14501193	10390.68299
23.13737991	10392.51999
23.13019719	10394.35791
23.12343166	10396.17051
23.11696447	10397.94558
23.11071537	10399.65842
23.1045814	10401.29831
23.09850078	10402.84888
23.09235201	10404.30986
23.08605499	10405.6721
23.07956724	10406.91654
23.07282876	10408.02727
23.06575404	10408.99828
23.05830642	10409.80715
23.05045417	10410.43684
23.04218362	10410.85022 
23.0334827	10411.06099
23.02439097	10411.06648
0 11500];


RPM = transpose(points(:,2));
TQ = transpose(points(:,1));

torque_fn = [RPM; TQ];
plot(RPM, TQ);