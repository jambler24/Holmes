from django.contrib import admin
from django.conf.urls import include, url
from django.conf.urls.static import static
from django.conf import settings

from holmes import views

urlpatterns = [
	url(r'^$', views.home, name='home'),
	url(r'^uploads/$', views.uploads, name='uploads'),
	url(r'^viewGene/(?P<gene_id>\w{0,50})/$', views.gene_view, name='geneView'),
	url(r'^searchGenes/$', views.gene_search, name='geneSearch'),
	url(r'^subNetworks/$', views.sub_graphs, name='subNetworks'),
	url(r'^varOverview/$', views.variant_overview, name='varOverview'),
	url(r'^subNetworks/Subnet_*?', views.sub_graph_detail, name='subNetworkDetail'),
	url(r'^coverageSummary/$', views.coverage_summary, name='coverageSummary'),
	url(r'^coverageSummaryGene/$', views.coverage_summary_gene, name='coverageSummaryGene'),
	url(r'^viewRegion/$', views.view_region, name='viewRegion'),
	url(r'^getGeneSubset/$', views.create_subset_file, name='getGeneSubset'),
	url('admin/', admin.site.urls),
] + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)


