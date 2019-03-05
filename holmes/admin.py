from django.contrib import admin

# Register your models here.

from .models import Experiment, CurrentSettings

admin.site.register(Experiment)
admin.site.register(CurrentSettings)

