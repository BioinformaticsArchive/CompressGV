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

	$('textarea')
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
			
			$('#results>table').not('#resultsTemplate').remove();
			
			var limits = [[4,novel-4],[novel,lines.length-5]];
			var ids = ['cross','novel'];
			for(var lim in limits){
				var table = $('#resultsTemplate').clone().attr('id', ids[lim]);
				var tbody = table.find('>tbody');
				for(var i=limits[lim][0]; i<limits[lim][1]; i++){
					var tr = $('<tr>');
					var vData = $.map(lines[i].split(/\s+/), function(v){
						return tr.append($('<td>').html(v));
					});
					tbody.append(tr);
				}
				if(lim==0){
					tbody.children().each(function(i,el){tbody.prepend(el)}); //reverse the order as the results are showing up reversed
				}
				$('#resultsTemplate').after(table.show());
			}
				
			$('#results>pre').html(lines[4]+'\n'+lines[novel]);
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
	
});
