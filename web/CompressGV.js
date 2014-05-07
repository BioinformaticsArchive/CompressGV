/*
Copyright 2013 Arran Schlosberg.

This file is part of https://github.com/aschlosberg/CompressGV (CompressGV)

    CompressGV is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CompressGV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CompressGV. If not, see <http://www.gnu.org/licenses/>.
*/

$(function(){
	$('#container>div').not(':first').hide();
	$('#nav>li:first').addClass('active');
	$('#year').html((new Date()).getFullYear());
	$('#resultsTemplate').hide();
	$('#runPHP').click(false);

	//will only work if prepare.sh has been run
	var ver = $('body').attr('data-version');
	if(ver!='master'){	
		$.getJSON('./uptodate.php', function(data){
			if(ver!=data.v){
				location.reload(true);
			}
		});
	}
	
	var toSection = function(section){
		var link = $('#nav>li>a[href=#'+section+']');
		if(!link){
			return false;
		}
		
		$('#nav>li.active').removeClass('active');
		link.parent().addClass('active');
	
		$('#container>div').hide();
		$('a[name='+section+']').parent().show();
	}
	
	var toHash = function(){
		var idx = location.href.indexOf("#");
		if(idx==-1){
			return false;
		}
		var hash = location.href.substring(idx+1);
		toSection(hash);
	}
	
	toHash();
	$(window).hashchange(toHash);

	$('textarea, input')
		.addClass('input-block-level')
		.attr('rows', 10)
		.attr('wrap', 'off')
		.css('font-family', 'Monaco,Menlo,Consolas,"Courier New",monospace');
		
	$('pre.example, #run>pre')
		.addClass('well')
		.addClass('pre-scrollable')
		.css('max-height', '80px');
		
	$('pre.example')
		.each(function(){
			var copy = $('<a href="#">Copy example to input</a>')
								.addClass('pull-right')
								.addClass('copy-example')
								.click(function(){
									$(this).parent().parent().find('textarea').val(
										$(this).parent().next('pre').html().replace(/&gt;/g, '>')
									);
									return false;
								});
								
			$(this).prev().append(copy);
		});
	
	var pdfSuffix = ['P', 'N'];
	var pdf = [];
	
	doPostTest = function(){
	
		if(pdf.length<2){
			for(var s in pdfSuffix){
				$.ajax({
					url : './pdf'+pdfSuffix[s]+'.json',
					async : false,
					success : function(data){
						pdf[s] = [];
						for(var d in data){
							pdf[s].push([-10 + d/100, data[d]]); // 2000 data points from -10 to 9.99 inclusive
						}
					}
				});
			}
		}
		
		var preTest = Math.max(0, Math.min(1, parseFloat($('#preTest').val())));
		var post = [];
		for(var i in pdf[0]){
			var inD = preTest * pdf[0][i][1];
			var inN = (1-preTest) * pdf[1][i][1];
			post[i] = [pdf[0][i][0], inD / (inD + inN)];
		}
	
		$.plot("#plot",
		[
			{ label: "Post-test Probability", color: 'black', data: post },
			{ label: "PDF (Deleterious)", color: 'red', data: pdf[0] },
			{ label: "PDF (Neutral)", color: 'blue', data: pdf[1] }
		],
		{
			xaxis : {
				min : -10,
				max : 10
			},
			yaxis : {
				min : 0,
				max : 1
			}
		});
		
		$('#probabilityTable tbody tr').each(function(){
			var GM = parseFloat($(this).find('td').eq(1).html());
			var idx = Math.floor((GM + 10) / 0.01);
			$(this).find('td').eq(2).html(Math.round(post[idx][1]*10000)/100+'%');
		});
	}

	doPostTest();
	$('#preTest').spinner({
		min : 0,
		max : 1,
		step : 0.01,
		change : doPostTest,
		spin : doPostTest
	});
		
	$('#runLink').click(function(){
		$('#run>div.alert').remove();
		var allOK = true;
		$('#run>pre').each(function(){
			var rel = $(this).attr('rel');
			var textarea = $('#'+rel+'>textarea');
			var val = $.trim(textarea.val());
			
			var lines = val.split(/\n/);
			for(var l in lines){
				lines[l] = $.trim(lines[l]).toUpperCase();
			}
			
			lines = lines.filter(function(n){return n}); //clear empty lines
			val = lines.join('\n');
			textarea.val(lines.join('\n'));
			
			var errors = [];
			var toMatch = rel=='msa'
								? /^([-ACDEFGHIKLMNPQRSTVWXY]+)$/
								: /^([-ACDEFGHIKLMNPQRSTVWXY][0-9]+[-ACDEFGHIKLMNPQRSTVWXY])$/;

			for(var l in lines){
				if(!(rel=='msa' && lines[l].substr(0,1)=='>') && !lines[l].match(toMatch)){
					errors.push(parseInt(l)+1);
				}
			}
			
			if((rel=='del' || rel=='neut') && lines.length<2){
				$(this).before('<div class="alert"><strong>WARNING!</strong> At least 2 nsSNPs required</div>');
				allOK = false;
			}
			else if(!lines.length){
				$(this).before('<div class="alert"><strong>WARNING!</strong> Empty</div>');
				allOK = false;
			}
			
			if(errors.length){
				$(this).before('<div class="alert"><strong>WARNING!</strong> Errors were found on lines: '+errors.join(', ') +'</div>');
				allOK = false;
			}
			
			$(this).html(lines.join('\n'));
		});
		
		if(allOK){
			$('#goForIt').removeClass('disabled');
		}
		else {
			$('#goForIt').addClass('disabled');
		}
	});

	$('#goForIt').click(function(){
		if($(this).hasClass('disabled')){
			return false;
		}
		
		$(this)
			.addClass('disabled')
			.html('Running...');

		var toPost = {};
		$('#run>pre').each(function(){
			var rel = $(this).attr('rel');
			var val = $(this).html();
			toPost[rel] = val;
		});

		$.post('./run.php', toPost, function(data){
			$('#goForIt')
				.removeClass('disabled')
				.html('Run');
				
			var lines = data.response.split('\n');
			var novel;
			for(var i=4; i<lines.length; i++){
				if(lines[i].substr(0,3)=='---'){
					novel = i+4;
					break;
				}
			}
			
			$('#results>*').not('#resultsTemplate').remove();
			var probTBody = $('#probabilityTable tbody');
			probTBody.find('>*').remove();
			
			var limits = [[4,novel-4],[novel,lines.length-5]];
			var ids = ['cross','novel'];
			var headings = ['Cross-validation', 'Novel nsSNPs'];
			for(var lim in limits){
				var table = $('#resultsTemplate').clone().attr('id', ids[lim]);
				$('#resultsTemplate').after(
					table.show()
				);
				table.before(
					$('<h3>').html(headings[lim])
				)

				var tbody = table.find('>tbody');
				for(var i=limits[lim][0]; i<limits[lim][1]; i++){
					var tr = $('<tr>');
					if(lim==1){
						var probTr = $('<tr>');
					}
					var vData = $.map(lines[i].split(/\s+/), function(v, i){
						var td = $('<td>').html(v);
						if(lim==1 && (i==0 || i==5)){
							probTr.append(td.clone());
						}
						return tr.append(td);
					});
					tbody.append(tr);
					if(lim==1){
						probTr.append($('<td>'));
						probTBody.append(probTr);
					}
				}

				if(lim==0){
					//calculate the AUC
					var gm = [];
					tbody.children().each(function(i,el){
						el = $(el);
						var cl;
						switch(el.children().eq(1).html()){
							case 'TN': case 'FP':
								cl = 'N';
								break;
							case 'TP': case 'FN':
								cl = 'P';
								break;
						}
						gm.push([cl, parseFloat(el.children().eq(5).html())]);
						tbody.prepend(el); //while we're at it... reverse the order as the results are showing up reversed
					});
					
					gm.sort(function(a, b){ return a[1]-b[1]; });
					var n = {N : 0, P : 0, product : 0};
					var beat = {N : 0, P : 0, sum : 0};
					for(var i=0; i<gm.length; i++){
						n[gm[i][0]]++;
						for(j=i+1; j<gm.length; j++){
							if(gm[i][0]!=gm[j][0]){
								beat[gm[i][0]]++;
							}
						}
					}
					
					n.product = n.N * n.P;
					beat.sum = beat.N + beat.P
					if(n.product==beat.sum){ //just checking :)
						table.before($('<h4>').html('ROC AUC: '+Math.round(beat.N / beat.sum*10000)/10000));
					}
				}
			}
				
			$('#results>pre').html(lines[4]+'\n'+lines[novel]);
			doPostTest();
			$('#resultLink').tab('show');
		}, 'json');

		return false;
	});
	
	$('#demo').click(function(){
		location.href = "#use";
		$('a.copy-example').click();
		$('#runLink').click();
		$('#goForIt').click();
		return false;
	});
	
	$('#compatible a').click(function(){
		var me = $(this).addClass('disabled');
		var rel = $(this).attr('rel');
		var tq = 'Thank you for your help, your response has been logged.';
		$.get($(this).attr('href'), function(){
			if(rel=='1'){
				alert(tq);
			}
			else {
				if(confirm(tq+' Would you like to lodge an issue and provide us with more information?')){
					location.href = 'https://github.com/aschlosberg/CompressGV/issues/new';
				}
			}
			me.removeClass('disabled');
		});
		return false;
	});
	
});
